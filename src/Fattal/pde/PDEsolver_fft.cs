using System;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using PdeSolver.Common;

namespace PdeSolver.FFT
{
    public static class PdeFftSolver
    {
        /// <summary>
        /// 고유벡터 공간의 데이터를 원래 공간으로 변환합니다. (T = EVy A EVx^tr)
        /// 입력 데이터 행렬 A를 변형시킵니다.
        /// </summary>
        public static void TransformEv2Normal(Array2Df A, Array2Df T)
        {
            int width = A.Width;
            int height = A.Height;

            // 올바른 변환을 위해 입력 값을 스케일링합니다.
            for (int y = 1; y < height - 1; y++)
            {
                for (int x = 1; x < width - 1; x++)
                {
                    A[x, y] *= 0.25f;
                }
            }

            for (int x = 1; x < width - 1; x++)
            {
                A[x, 0] *= 0.5f;
                A[x, height - 1] *= 0.5f;
            }

            for (int y = 1; y < height - 1; y++)
            {
                A[0, y] *= 0.5f;
                A[width - 1, y] *= 0.5f;
            }

            // 2D 이산 코사인 변환(DCT) 실행
            ExecuteFft2D(A, T, height, width);
        }

        /// <summary>
        /// 원래 공간의 데이터를 고유벡터 공간으로 변환합니다. (T = EVy^-1 * A * (EVx^-1)^tr)
        /// </summary>
        public static void TransformNormal2Ev(Array2Df A, Array2Df T)
        {
            int width = A.Width;
            int height = A.Height;

            // 2D 이산 코사인 변환(DCT) 실행
            ExecuteFft2D(A, T, height, width);

            // 올바른 변환을 위해 출력 행렬을 스케일링합니다.
            float scaleFactor = 1.0f / ((height - 1) * (width - 1));
            
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    T[x, y] *= scaleFactor;
                }
            }

            for (int x = 0; x < width; x++)
            {
                T[x, 0] *= 0.5f;
                T[x, height - 1] *= 0.5f;
            }

            for (int y = 0; y < height; y++)
            {
                T[0, y] *= 0.5f;
                T[width - 1, y] *= 0.5f;
            }
        }

        /// <summary>
        /// 1차원 라플라스 연산자의 고유값을 반환합니다.
        /// </summary>
        private static double[] GetLambda(int n)
        {
            double[] v = new double[n];
            for (int i = 0; i < n; i++)
            {
                double sinVal = Math.Sin((double)i / (2 * (n - 1)) * Math.PI);
                v[i] = -4.0 * sinVal * sinVal;
            }
            return v;
        }

        /// <summary>
        /// 편미분 방정식의 해가 존재하도록 경계 조건을 호환되게 조정합니다.
        /// </summary>
        public static void MakeCompatibleBoundary(Array2Df F)
        {
            int width = F.Width;
            int height = F.Height;

            double sum = 0.0;
            for (int y = 1; y < height - 1; y++)
            {
                for (int x = 1; x < width - 1; x++) sum += F[x, y];
            }

            for (int x = 1; x < width - 1; x++)
            {
                sum += 0.5 * (F[x, 0] + F[x, height - 1]);
            }

            for (int y = 1; y < height - 1; y++)
            {
                sum += 0.5 * (F[0, y] + F[width - 1, y]);
            }

            sum += 0.25 * (F[0, 0] + F[0, height - 1] + F[width - 1, 0] + F[width - 1, height - 1]);

            float add = (float)(-sum / (height + width - 3));

            for (int x = 0; x < width; x++)
            {
                F[x, 0] += add;
                F[x, height - 1] += add;
            }

            for (int y = 1; y < height - 1; y++)
            {
                F[0, y] += add;
                F[width - 1, y] += add;
            }
        }

        /// <summary>
        /// 노이만 경계 조건을 사용하여 라플라스 방정식 $U = F$ 를 풉니다.
        /// adjust_bound가 true이면 솔루션이 존재하도록 경계값을 조정합니다.
        /// </summary>
        public static void SolvePdeFft(Array2Df F, Array2Df U, Array2Df F_tr, bool adjustBound, Action<int> progressCallback = null)
        {
            progressCallback?.Invoke(20);

            int width = F.Width;
            int height = F.Height;

            if (adjustBound)
            {
                MakeCompatibleBoundary(F);
            }

            // F를 고유벡터 공간으로 변환
            TransformNormal2Ev(F, F_tr);
            
            progressCallback?.Invoke(50);

            double[] l1 = GetLambda(height);
            double[] l2 = GetLambda(width);

            // 고유벡터 공간에서 해를 계산 (병렬 처리 적용)
            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    F_tr[x, y] = (float)(F_tr[x, y] / (l1[y] + l2[x]));
                }
            });

            // 원점의 값은 0으로 설정 (해에 상수를 더하는 것과 동일함)
            F_tr[0, 0] = 0.0f; 

            progressCallback?.Invoke(55);

            // 해를 원래 공간으로 역변환
            TransformEv2Normal(F_tr, U);
            
            progressCallback?.Invoke(85);

            // 양수 값을 제거하기 위해 상수를 제거 (수치적 안정성을 위해)
            float max = float.MinValue;
            int length = U.Length;
            
            for (int i = 0; i < length; i++)
            {
                float val = U.Get(i);
                if (max < val)
                {
                    max = val;
                }
            }

            for (int i = 0; i < length; i++)
            {
                U.Set(i, U.Get(i) - max);
            }

            progressCallback?.Invoke(90);
        }

        /// <summary>
        /// 내부 점에 대한 (Laplace U - F)의 잔차(Residual) 노름(Norm)을 반환합니다.
        /// 솔버의 정확도를 확인할 때 사용됩니다.
        /// </summary>
        public static float ResidualPde(Array2Df U, Array2Df F)
        {
            int width = U.Width;
            int height = U.Height;

            double res = 0.0;
            for (int y = 1; y < height - 1; y++)
            {
                for (int x = 1; x < width - 1; x++)
                {
                    double laplace = -4.0 * U[x, y] + U[x - 1, y] + U[x + 1, y] + U[x, y - 1] + U[x, y + 1];
                    double diff = laplace - F[x, y];
                    res += diff * diff;
                }
            }
            return (float)Math.Sqrt(res);
        }

        /// <summary>
        /// 외부 FFT 라이브러리 연동을 위한 인터페이스 지점입니다.
        /// C# 프로젝트 환경에 맞춰 Fftw.Net 또는 MathNet.Numerics 등의 DCT-1 (REDFT00) 기능을 연결해야 합니다.
        /// </summary>
        private static void ExecuteFft2D(Array2Df input, Array2Df output, int height, int width)
        {
            // TODO: 이곳에 외부 FFT 라이브러리 호출 코드를 작성하십시오.
            // 예시 (가상 코드):
            // FftwWrapper.ExecuteR2R(height, width, input.Data, output.Data, FftwTransformType.REDFT00);
            throw new NotImplementedException("FFT 라이브러리 구현체가 필요합니다.");
        }
    }
}