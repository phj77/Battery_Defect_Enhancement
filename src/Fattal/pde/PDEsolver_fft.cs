using System;
using System.Threading.Tasks;
using SharpFFTW.Single;
using SharpFFTW;

namespace PdeSolver.FFT
{
    /// <summary>
    /// Direct Poisson solver using the discrete cosine transform.
    /// </summary>
    public static class PdeFftSolver
    {
        /// <summary>
        /// 1D 라플라스 연산자의 고유값(Eigenvalues)을 반환합니다.
        /// </summary>
        private static float[] GetLambda(int n)
        {
            float[] v = new float[n];
            for (int i = 0; i < n; i++)
            {
                v[i] = (float)(-4.0 * Math.Pow(Math.Sin((double)i / (2.0 * (n - 1)) * Math.PI), 2));
            }
            return v;
        }

        /// <summary>
        /// 해(Solution)가 존재할 수 있도록 경계 조건을 호환되게 조정합니다.
        /// </summary>
        private static void MakeCompatibleBoundary(float[] F, int width, int height)
        {
            double sum = 0.0;

            // 내부 영역 합산
            for (int y = 1; y < height - 1; y++)
            {
                int rowOffset = y * width;
                for (int x = 1; x < width - 1; x++)
                {
                    sum += F[rowOffset + x];
                }
            }

            // 가장자리 합산
            for (int x = 1; x < width - 1; x++)
            {
                sum += 0.5 * (F[x] + F[(height - 1) * width + x]);
            }
            for (int y = 1; y < height - 1; y++)
            {
                int rowOffset = y * width;
                sum += 0.5 * (F[rowOffset] + F[rowOffset + (width - 1)]);
            }

            // 모서리 합산
            sum += 0.25 * (F[0] + F[(height - 1) * width] + F[width - 1] + F[(height - 1) * width + (width - 1)]);

            double add = -sum / (height + width - 3);
            float fAdd = (float)add;

            // 경계값 조정
            for (int x = 0; x < width; x++)
            {
                F[x] += fAdd;
                F[(height - 1) * width + x] += fAdd;
            }
            for (int y = 1; y < height - 1; y++)
            {
                int rowOffset = y * width;
                F[rowOffset] += fAdd;
                F[rowOffset + (width - 1)] += fAdd;
            }
        }

        /// <summary>
        /// 일반 공간에서 고유벡터 공간으로 변환합니다. (T = EVy^-1 * A * (EVx^-1)^tr)
        /// </summary>
        private static void TransformNormal2Ev(float[] A, float[] T, int width, int height)
        {
            // 1. 스레드 환경 초기화 (성공 시 0이 아닌 값 반환)
            if (NativeMethods.fftwf_init_threads() == 0)
            {
                throw new Exception("fail to FFTW initialization!");
            }

            // 2. 사용할 CPU 스레드 수 지정 (12~ is optimal for local pc)
            NativeMethods.fftwf_plan_with_nthreads(14);

            Console.WriteLine($"fft plan start=================: {GlobalTimer.ElapsedSeconds:F2}s");
            using (var inArray = new RealArray(A))
            using (var outArray = new RealArray(T.Length))
            using (var plan = Plan.Create2(height, width, inArray, outArray, Transform.REDFT00, Transform.REDFT00, Options.Estimate))
            {
                plan.Execute();
                outArray.CopyTo(T);
            }
            Console.WriteLine($"fft plan over=================: {GlobalTimer.ElapsedSeconds:F2}s");

            float scale = 1.0f / ((height - 1) * (width - 1));

            // 스케일링 병렬 처리
            Parallel.For(0, height, y =>
            {
                int rowOffset = y * width;
                for (int x = 0; x < width; x++)
                {
                    T[rowOffset + x] *= scale;
                }
            });

            // 가장자리 추가 스케일링
            for (int x = 0; x < width; x++)
            {
                T[x] *= 0.5f;
                T[(height - 1) * width + x] *= 0.5f;
            }
            for (int y = 0; y < height; y++)
            {
                int rowOffset = y * width;
                T[rowOffset] *= 0.5f;
                T[rowOffset + (width - 1)] *= 0.5f;
            }
        }

        /// <summary>
        /// 고유벡터 공간에서 일반 공간으로 변환합니다. (T = EVy A EVx^tr)
        /// 주의: 입력 배열 A의 데이터가 변형됩니다.
        /// </summary>
        private static void TransformEv2Normal(float[] A, float[] T, int width, int height)
        {
            // 입력 데이터 스케일링 병렬 처리
            Parallel.For(1, height - 1, y =>
            {
                int rowOffset = y * width;
                for (int x = 1; x < width - 1; x++)
                {
                    A[rowOffset + x] *= 0.25f;
                }
            });

            for (int x = 1; x < width - 1; x++)
            {
                A[x] *= 0.5f;
                A[(height - 1) * width + x] *= 0.5f;
            }
            for (int y = 1; y < height - 1; y++)
            {
                int rowOffset = y * width;
                A[rowOffset] *= 0.5f;
                A[rowOffset + (width - 1)] *= 0.5f;
            }

            Console.WriteLine($"fft plan start=================: {GlobalTimer.ElapsedSeconds:F2}s");
            using (var inArray = new RealArray(A))
            using (var outArray = new RealArray(T.Length))
            using (var plan = Plan.Create2(height, width, inArray, outArray, Transform.REDFT00, Transform.REDFT00, Options.Estimate))
            {
                plan.Execute();
                outArray.CopyTo(T);
            }
            Console.WriteLine($"fft plan over=================: {GlobalTimer.ElapsedSeconds:F2}s");
        }

        /// <summary>
        /// 노이만 경계 조건을 가진 라플라스 방정식 U = F 를 풉니다.
        /// </summary>
        public static void SolvePdeFft(float[] F, float[] U, float[] F_tr, int width, int height, bool adjustBound = true, float hpfSigma = 0.0f)
        {
            if (F.Length != width * height || U.Length != width * height || F_tr.Length != width * height)
            {
                throw new ArgumentException("Array lengths must match width * height.");
            }

            if (adjustBound)
            {
                Console.WriteLine($"Making compatible boundary start: {GlobalTimer.ElapsedSeconds:F2}s");
                MakeCompatibleBoundary(F, width, height);
            }

            Console.WriteLine($"Transform Divergence to eigenvector space start: {GlobalTimer.ElapsedSeconds:F2}s");
            // 1. F를 고유벡터 공간으로 변환
            TransformNormal2Ev(F, F_tr, width, height);

            float[] l1 = GetLambda(height);
            float[] l2 = GetLambda(width);

            Console.WriteLine($"calculate in eigenvector space start: {GlobalTimer.ElapsedSeconds:F2}s");
            // 2. 고유벡터 공간에서의 계산 및 High-Pass Filter 적용
            if (hpfSigma > 0.0f)
            {
                double sigma2 = 2.0 * hpfSigma * hpfSigma;

                // 속도 최적화: x축에 대한 제곱 값을 사전에 계산하여 내부 루프 연산량 감소
                double[] kx2 = new double[width];
                for (int x = 0; x < width; x++)
                {
                    double kx = (double)x / width;
                    kx2[x] = kx * kx;
                }

                Parallel.For(0, height, y =>
                {
                    int rowOffset = y * width;
                    float ly = l1[y];

                    // y축 계산값 캐싱
                    double ky = (double)y / height;
                    double ky2 = ky * ky;

                    for (int x = 0; x < width; x++)
                    {
                        if (x == 0 && y == 0) continue;

                        double hFilter = 1.0 - Math.Exp(-(ky2 + kx2[x]) / sigma2);
                        F_tr[rowOffset + x] = (float)((F_tr[rowOffset + x] / (ly + l2[x])) * hFilter);
                    }
                });
            }
            else
            {
                // 기존 계산 (가장 부하가 큰 연산이므로 Parallel 처리)
                Parallel.For(0, height, y =>
                {
                    int rowOffset = y * width;
                    float ly = l1[y];
                    for (int x = 0; x < width; x++)
                    {
                        // (0,0)에서의 0 나누기 방지
                        if (x == 0 && y == 0) continue;
                        F_tr[rowOffset + x] /= (ly + l2[x]);
                    }
                });
            }

            // 상수를 상쇄하기 위해 설정
            F_tr[0] = 0f;

            Console.WriteLine($"Inverse transform to time domain start: {GlobalTimer.ElapsedSeconds:F2}s");
            // 3. F_tr을 일반 공간으로 변환하여 U 도출
            TransformEv2Normal(F_tr, U, width, height);

            Console.WriteLine($"subtract maximum value: {GlobalTimer.ElapsedSeconds:F2}s");
            // 4. 결과값에서 최대값을 찾아 빼줌 (로그 공간 연산의 수치적 안정성을 위함)
            float max = float.MinValue;
            for (int i = 0; i < U.Length; i++)
            {
                if (U[i] > max)
                {
                    max = U[i];
                }
            }

            Parallel.For(0, U.Length, i =>
            {
                U[i] -= max;
            });
        }

        /// <summary>
        /// 정확도 검증을 위한 내부 점들의 (Laplace U - F) 노름(norm) 반환 함수
        /// </summary>
        public static float ResidualPde(float[] U, float[] F, int width, int height)
        {
            double res = 0.0;
            object lockObj = new object();

            // Thread-Local 변수를 활용한 고속 Map-Reduce 병렬 처리
            Parallel.For(1, height - 1, () => 0.0, (y, loopState, localRes) =>
            {
                int rowOffset = y * width;
                int upOffset = (y - 1) * width;
                int downOffset = (y + 1) * width;

                for (int x = 1; x < width - 1; x++)
                {
                    double laplace = -4.0 * U[rowOffset + x] +
                                     U[rowOffset + x - 1] + U[rowOffset + x + 1] +
                                     U[upOffset + x] + U[downOffset + x];

                    double diff = laplace - F[rowOffset + x];
                    localRes += diff * diff;
                }
                return localRes;
            },
            localRes =>
            {
                lock (lockObj) res += localRes;
            });

            return (float)Math.Sqrt(res);
        }
    }
}