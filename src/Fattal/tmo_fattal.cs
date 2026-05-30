using System;
using System.Threading.Tasks;
using PdeSolver.FFT;
using PdeSolver.Multigrid;
using PdeSolver.Common;

namespace FattalToneMapping
{
    public static class TmoFattal02
    {
        private static void DownSample(Array2Df A, Array2Df B)
        {
            int width = B.Cols;
            int height = B.Rows;

            Parallel.For(0, height, y =>
            {
                int y2 = y * 2;
                for (int x = 0; x < width; x++)
                {
                    int x2 = x * 2;
                    float p = A[x2, y2] + A[x2 + 1, y2] + A[x2, y2 + 1] + A[x2 + 1, y2 + 1];
                    B[x, y] = p * 0.25f;
                }
            });
        }

        private static void GaussianBlur(Array2Df I, Array2Df L)
        {
            int width = I.Cols;
            int height = I.Rows;

            if (width < 3 || height < 3)
            {
                if (!ReferenceEquals(I, L))
                {
                    Array.Copy(I.Data, L.Data, I.Length);
                }
                return;
            }

            Array2Df T = new Array2Df(width, height);

            // --- X축 블러 ---
            Parallel.For(0, height, y =>
            {
                for (int x = 1; x < width - 1; x++)
                {
                    float t = 2f * I[x, y] + I[x - 1, y] + I[x + 1, y];
                    T[x, y] = t * 0.25f;
                }
                T[0, y] = (3f * I[0, y] + I[1, y]) * 0.25f;
                T[width - 1, y] = (3f * I[width - 1, y] + I[width - 2, y]) * 0.25f;
            });

            // --- Y축 블러 ---
            Parallel.For(0, width, x =>
            {
                for (int y = 1; y < height - 1; y++)
                {
                    float t = 2f * T[x, y] + T[x, y - 1] + T[x, y + 1];
                    L[x, y] = t * 0.25f;
                }
                L[x, 0] = (3f * T[x, 0] + T[x, 1]) * 0.25f;
                L[x, height - 1] = (3f * T[x, height - 1] + T[x, height - 2]) * 0.25f;
            });
        }

        private static void CreateGaussianPyramids(Array2Df H, Array2Df[] pyramids, int nlevels)
        {
            int width = H.Cols;
            int height = H.Rows;

            Array2Df L = new Array2Df(width, height);
            GaussianBlur(pyramids[0], L);

            for (int k = 1; k < nlevels; k++)
            {
                width /= 2;
                height /= 2;
                pyramids[k] = new Array2Df(width, height);
                DownSample(L, pyramids[k]);

                if (k < nlevels - 1)
                {
                    L = new Array2Df(width, height); // 이전 L 참조 버리기
                    GaussianBlur(pyramids[k], L);
                }
            }
        }

        private static float CalculateGradients(Array2Df H, Array2Df G, int k)
        {
            int width = H.Cols;
            int height = H.Rows;
            float divider = MathF.Pow(2.0f, k + 1);

            object lockObj = new object();
            double avgGrad = 0.0;

            Parallel.For(0, height, () => 0.0, (y, loopState, localSum) =>
            {
                for (int x = 0; x < width; x++)
                {
                    int w = (x == 0 ? 0 : x - 1);
                    int n = (y == 0 ? 0 : y - 1);
                    int s = (y + 1 == height ? y : y + 1);
                    int e = (x + 1 == width ? x : x + 1);

                    float gx = (H[w, y] - H[e, y]) / divider;
                    float gy = (H[x, s] - H[x, n]) / divider;

                    float mag = MathF.Sqrt(gx * gx + gy * gy);
                    G[x, y] = mag;
                    localSum += mag;
                }
                return localSum;
            }, localSum => { lock (lockObj) avgGrad += localSum; });

            return (float)(avgGrad / (width * height));
        }

        private static void UpSample(Array2Df A, Array2Df B)
        {
            int width = B.Cols;
            int height = B.Rows;
            int awidth = A.Cols;
            int aheight = A.Rows;

            Parallel.For(0, height, y =>
            {
                int ay = Math.Min((int)(y * 0.5f), aheight - 1);
                for (int x = 0; x < width; x++)
                {
                    int ax = Math.Min((int)(x * 0.5f), awidth - 1);
                    B[x, y] = A[ax, ay];
                }
            });
        }

        private static void CalculateFiMatrix(Array2Df FI, Array2Df[] gradients, float[] avgGrad,
                                              int nlevels, int detailLevel, float alfa, float beta,
                                              float noise, bool newFattal)
        {
            Array2Df[] fi = new Array2Df[nlevels];

            int width = gradients[nlevels - 1].Cols;
            int height = gradients[nlevels - 1].Rows;
            fi[nlevels - 1] = new Array2Df(width, height);

            if (newFattal)
            {
                for (int k = 0; k < width * height; k++) fi[nlevels - 1].Data[k] = 1.0f;
            }

            for (int k = nlevels - 1; k >= 0; k--)
            {
                width = gradients[k].Cols;
                height = gradients[k].Rows;

                if (k >= detailLevel || k == nlevels - 1 || !newFattal)
                {
                    Parallel.For(0, height, y =>
                    {
                        for (int x = 0; x < width; x++)
                        {
                            float grad = MathF.Max(gradients[k][x, y], 1e-4f);
                            float a = alfa * avgGrad[k];
                            float value = MathF.Pow((grad + noise) / a, beta - 1.0f);

                            if (newFattal)
                                fi[k][x, y] *= value;
                            else
                                fi[k][x, y] = value;
                        }
                    });
                }

                if (k > 1)
                {
                    width = gradients[k - 1].Cols;
                    height = gradients[k - 1].Rows;
                    fi[k - 1] = new Array2Df(width, height);
                }
                else
                {
                    fi[0] = FI;
                }

                if (k > 0 && newFattal)
                {
                    UpSample(fi[k], fi[k - 1]);
                    GaussianBlur(fi[k - 1], fi[k - 1]);
                }
            }
        }

        private static void FindMinMaxPercentile(float[] data, float cutMin, out float minVal, float cutMax, out float maxVal)
        {
            float[] sorted = (float[])data.Clone();
            Array.Sort(sorted);

            int minIdx = (int)(sorted.Length * cutMin);
            int maxIdx = (int)(sorted.Length * cutMax);

            minIdx = Math.Clamp(minIdx, 0, sorted.Length - 1);
            maxIdx = Math.Clamp(maxIdx, 0, sorted.Length - 1);

            minVal = sorted[minIdx];
            maxVal = sorted[maxIdx];
        }

        public static void Process(int width, int height, Array2Df Y, Array2Df L,
                                   float alfa, float beta, float noise, bool newFattal,
                                   bool fftSolver, int detailLevel, Action<int> progressCallback)
        {
            const float blackPoint = 0.1f;
            const float whitePoint = 0.5f;
            const float gamma = 1.0f;

            detailLevel = Math.Clamp(detailLevel, 0, 3);
            progressCallback?.Invoke(2);

            int msize = fftSolver ? 8 : 32;
            int size = width * height;

            // 1. 최대/최소 루미넌스 탐색
            float minLum = Y.Data[0];
            float maxLum = Y.Data[0];

            for (int i = 0; i < size; i++)
            {
                if (Y.Data[i] < minLum) minLum = Y.Data[i];
                if (Y.Data[i] > maxLum) maxLum = Y.Data[i];
            }

            // 2. 로그 도메인으로 변환 (H)
            Array2Df H = new Array2Df(width, height);
            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    H[x, y] = MathF.Log(100f * Y[x, y] / maxLum + 1e-4f);
                }
            });

            progressCallback?.Invoke(4);

            // 3. 가우시안 피라미드 생성
            int mins = Math.Min(width, height);
            int nlevels = 0;
            while (mins >= msize)
            {
                nlevels++;
                mins /= 2;
            }
            if (nlevels == 0) nlevels = 1;

            Array2Df[] pyramids = new Array2Df[nlevels];
            pyramids[0] = H;
            CreateGaussianPyramids(H, pyramids, nlevels);

            progressCallback?.Invoke(8);

            // 4. 그래디언트 피라미드 계산
            Array2Df[] gradients = new Array2Df[nlevels];
            float[] avgGrad = new float[nlevels];
            for (int k = 0; k < nlevels; k++)
            {
                gradients[k] = new Array2Df(pyramids[k].Cols, pyramids[k].Rows);
                avgGrad[k] = CalculateGradients(pyramids[k], gradients[k], k);
            }

            progressCallback?.Invoke(12);

            // 5. 감쇠 매트릭스(FI) 계산
            Array2Df FI = new Array2Df(width, height);
            CalculateFiMatrix(FI, gradients, avgGrad, nlevels, detailLevel, alfa, beta, noise, newFattal);

            progressCallback?.Invoke(16);

            // 6. 감쇠된 그래디언트 필드 생성
            Array2Df Gx = new Array2Df(width, height);
            Array2Df Gy = new Array2Df(width, height);

            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    if (fftSolver)
                    {
                        int yp1 = y + 1 >= height ? height - 2 : y + 1;
                        int xp1 = x + 1 >= width ? width - 2 : x + 1;
                        Gx[x, y] = (H[xp1, y] - H[x, y]) * 0.5f * (FI[xp1, y] + FI[x, y]);
                        Gy[x, y] = (H[x, yp1] - H[x, y]) * 0.5f * (FI[x, yp1] + FI[x, y]);
                    }
                    else
                    {
                        int s = y + 1 == height ? y : y + 1;
                        int e = x + 1 == width ? x : x + 1;
                        Gx[x, y] = (H[e, y] - H[x, y]) * FI[x, y];
                        Gy[x, y] = (H[x, s] - H[x, y]) * FI[x, y];
                    }
                }
            });

            progressCallback?.Invoke(18);

            // 7. Divergence 계산
            Array2Df DivG = new Array2Df(width, height);
            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    DivG[x, y] = Gx[x, y] + Gy[x, y];
                    if (x > 0) DivG[x, y] -= Gx[x - 1, y];
                    if (y > 0) DivG[x, y] -= Gy[x, y - 1];

                    if (fftSolver)
                    {
                        if (x == 0) DivG[x, y] += Gx[x, y];
                        if (y == 0) DivG[x, y] += Gy[x, y];
                    }
                }
            });

            progressCallback?.Invoke(20);

            // 8. 포아송 방정식 풀이 (PDE Solver)
            Array2Df U = new Array2Df(width, height);
            if (fftSolver)
            {
                // 수정된 부분: Array2Df 인스턴스 대신 내부 Data(float[])와 width, height를 전달합니다.
                // PdeFftSolver.SolvePdeFft 메서드 서명에 맞추어 인자를 구성했습니다.
                PdeFftSolver.SolvePdeFft(DivG.Data, U.Data, Gx.Data, width, height, false, 0.007f);
            }
            else
            {
                PdeMultigridSolver.SolvePdeMultigrid(DivG, U, progressCallback);
            }

            progressCallback?.Invoke(90);

            // 9. 압축된 이미지 복원 (지수 함수)
            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    L[x, y] = MathF.Exp(gamma * U[x, y]);
                }
            });

            progressCallback?.Invoke(95);

            // 10. 백분위수를 이용한 정규화 및 클리핑
            float cutMin = 0.01f * blackPoint;
            float cutMax = 1.0f - 0.01f * whitePoint;

            FindMinMaxPercentile(L.Data, cutMin, out minLum, cutMax, out maxLum);

            Parallel.For(0, size, i =>
            {
                float val = (L.Data[i] - minLum) / (maxLum - minLum);
                L.Data[i] = val <= 0.0f ? 0.0f : val;
            });

            progressCallback?.Invoke(96);
        }
    }
}