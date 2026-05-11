using System;
using System.Buffers;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using PdeSolver.Common;

namespace PdeSolver.Multigrid
{
    public static class PdeMultigridSolver
    {
        // Multigrid 솔버 튜닝 파라미터
        private const int MODYF = 0;
        private const int MINS = 16;
        private const int SMOOTH_IT = 1;
        private const int BCG_STEPS = 20;
        private const float BCG_TOL = 1e-3f;
        private const int V_CYCLE = 2;
        private const bool BCG_POST_IMPROVE = false;
        private const int BCG_POST_STEPS = 1000;
        private const float BCG_POST_TOL = 1e-7f;
        private const int OMP_THRESHOLD = 1_000_000;

        /// <summary>
        /// 고해상도 그리드의 데이터를 저해상도 그리드로 축소(Restriction)합니다.
        /// </summary>
        private static void Restrict(Array2Df input, Array2Df output)
        {
            float inRows = input.Rows;
            float inCols = input.Cols;
            int outRows = output.Rows;
            int outCols = output.Cols;

            float dx = inCols / outCols;
            float dy = inRows / outRows;
            const float filterSize = 0.5f;

            float sy = dy / 2f - 0.5f;
            for (int y = 0; y < outRows; y++, sy += dy)
            {
                float sx = dx / 2f - 0.5f;
                for (int x = 0; x < outCols; x++, sx += dx)
                {
                    float pixVal = 0;
                    float w = 0;

                    float startX = MathF.Max(0, MathF.Ceiling(sx - dx * filterSize));
                    float endX = MathF.Min(MathF.Floor(sx + dx * filterSize), inCols - 1);
                    float startY = MathF.Max(0, MathF.Ceiling(sy - dy * filterSize));
                    float endY = MathF.Min(MathF.Floor(sy + dy * filterSize), inRows - 1);

                    for (float ix = startX; ix <= endX; ix++)
                    {
                        for (float iy = startY; iy <= endY; iy++)
                        {
                            pixVal += input[(int)ix, (int)iy];
                            w += 1f;
                        }
                    }
                    output[x, y] = pixVal / w;
                }
            }
        }

        /// <summary>
        /// 저해상도 그리드의 데이터를 고해상도 그리드로 보간(Prolongation)합니다.
        /// </summary>
        private static void Prolongate(Array2Df input, Array2Df output)
        {
            float dx = (float)input.Cols / output.Cols;
            float dy = (float)input.Rows / output.Rows;

            int outRows = output.Rows;
            int outCols = output.Cols;
            float inRows = input.Rows;
            float inCols = input.Cols;
            const float filterSize = 1f;

            float sy = -dy / 2f;
            for (int y = 0; y < outRows; y++, sy += dy)
            {
                float sx = -dx / 2f;
                for (int x = 0; x < outCols; x++, sx += dx)
                {
                    float pixVal = 0;
                    float weight = 0;

                    float startX = MathF.Max(0, MathF.Ceiling(sx - filterSize));
                    float endX = MathF.Min(MathF.Floor(sx + filterSize), inCols - 1);
                    float startY = MathF.Max(0, MathF.Ceiling(sy - filterSize));
                    float endY = MathF.Min(MathF.Floor(sy + filterSize), inRows - 1);

                    for (float ix = startX; ix <= endX; ix++)
                    {
                        for (float iy = startY; iy <= endY; iy++)
                        {
                            float fx = MathF.Abs(sx - ix);
                            float fy = MathF.Abs(sy - iy);
                            float fval = (1f - fx) * (1f - fy);

                            pixVal += input[(int)ix, (int)iy] * fval;
                            weight += fval;
                        }
                    }
                    output[x, y] = pixVal / weight;
                }
            }
        }

        /// <summary>
        /// 최하위(Coarsest) 그리드에서의 해를 구합니다.
        /// </summary>
        private static void ExactSolution(Array2Df F, Array2Df U)
        {
            U.Reset();
        }

        /// <summary>
        /// Biconjugate Gradient Method를 이용해 현재 해(U)를 평활화(Smoothing)합니다.
        /// </summary>
        private static void Smooth(Array2Df U, Array2Df F)
        {
            int n = U.Length;
            LinBcg(n, F.Data, U.Data, BCG_TOL, BCG_STEPS, out _, out _, U.Rows, U.Cols);
        }

        /// <summary>
        /// 편미분 방정식의 잔여 오차(Defect/Residual)를 계산합니다. D = F - L*U
        /// </summary>
        private static void CalculateDefect(Array2Df D, Array2Df U, Array2Df F)
        {
            int sx = F.Cols;
            int sy = F.Rows;

            for (int y = 0; y < sy; y++)
            {
                for (int x = 0; x < sx; x++)
                {
                    int w = (x == 0 ? 0 : x - 1);
                    int n = (y == 0 ? 0 : y - 1);
                    int s = (y + 1 == sy ? y : y + 1);
                    int e = (x + 1 == sx ? x : x + 1);

                    D[x, y] = F[x, y] - (U[e, y] + U[w, y] + U[x, n] + U[x, s] - 4.0f * U[x, y]);
                }
            }
        }

        /// <summary>
        /// 현재 해에 보정값(Correction)을 더합니다.
        /// </summary>
        private static void AddCorrection(Array2Df U, Array2Df C)
        {
            int n = U.Length;
            VAdds(U.Data, 1.0f, C.Data, U.Data, n);
        }

        /// <summary>
        /// 두 Array2Df 간의 단순 메모리 복사를 수행합니다.
        /// </summary>
        private static void Copy(Array2Df src, Array2Df dst)
        {
            Array.Copy(src.Data, dst.Data, src.Length);
        }

        /// <summary>
        /// Full Multigrid 알고리즘의 메인 루틴입니다.
        /// 다중 해상도 그리드를 구성하고, V-Cycle을 반복하며 PDE의 해를 구합니다.
        /// </summary>
        public static void SolvePdeMultigrid(Array2Df F, Array2Df U, Action<int> progressCallback)
        {
            int xmax = F.Cols;
            int ymax = F.Rows;
            int levels = 0;
            int mins = Math.Min(xmax, ymax);

            while (mins >= MINS)
            {
                levels++;
                mins = mins / 2 + MODYF;
            }

            Array2Df[] RHS = new Array2Df[levels + 1];
            Array2Df[] IU = new Array2Df[levels + 1];
            Array2Df[] VF = new Array2Df[levels + 1];

            VF[0] = new Array2Df(xmax, ymax);
            RHS[0] = F;
            IU[0] = new Array2Df(xmax, ymax);
            Copy(U, IU[0]);

            int sx = xmax;
            int sy = ymax;

            for (int k = 0; k < levels; k++)
            {
                sx = sx / 2 + MODYF;
                sy = sy / 2 + MODYF;

                RHS[k + 1] = new Array2Df(sx, sy);
                IU[k + 1] = new Array2Df(sx, sy);
                VF[k + 1] = new Array2Df(sx, sy);

                Restrict(RHS[k], RHS[k + 1]);
            }

            ExactSolution(RHS[levels], IU[levels]);

            for (int k = levels - 1; k >= 0; k--)
            {
                progressCallback?.Invoke(20 + 70 * (levels - k) / (levels + 1));
                Prolongate(IU[k + 1], IU[k]);
                Copy(RHS[k], VF[k]);

                for (int cycle = 0; cycle < V_CYCLE; cycle++)
                {
                    for (int k2 = k; k2 < levels; k2++)
                    {
                        if (k2 != k) IU[k2].Reset();

                        for (int i = 0; i < SMOOTH_IT; i++) Smooth(IU[k2], VF[k2]);

                        var D = new Array2Df(IU[k2].Cols, IU[k2].Rows);
                        CalculateDefect(D, IU[k2], VF[k2]);
                        Restrict(D, VF[k2 + 1]);
                    }

                    ExactSolution(VF[levels], IU[levels]);

                    for (int k2 = levels - 1; k2 >= k; k2--)
                    {
                        var C = new Array2Df(IU[k2].Cols, IU[k2].Rows);
                        Prolongate(IU[k2 + 1], C);
                        AddCorrection(IU[k2], C);

                        for (int i = 0; i < SMOOTH_IT; i++) Smooth(IU[k2], VF[k2]);
                    }
                }
            }

            Copy(IU[0], U);

            if (BCG_POST_IMPROVE)
            {
                LinBcg(xmax * ymax, F.Data, U.Data, BCG_POST_TOL, BCG_POST_STEPS, out _, out _, ymax, xmax);
            }

            progressCallback?.Invoke(90);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ASolve(float[] b, float[] x, int n)
        {
            int i = 0;
            int vectorLength = Vector<float>.Count;
            if (n >= vectorLength)
            {
                var vNeg4 = new Vector<float>(-4.0f);
                for (; i <= n - vectorLength; i += vectorLength)
                {
                    (new Vector<float>(b, i) * vNeg4).CopyTo(x, i);
                }
            }
            for (; i < n; i++) x[i] = -4.0f * b[i];
        }

        /// <summary>
        /// 2D 라플라스 연산자와 벡터 x의 행렬-벡터 곱을 수행합니다. (Neumann 경계 조건 포함)
        /// </summary>
        private static void ATimes(float[] x, float[] res, int rows, int cols)
        {
            // 중앙 영역
            if (rows * cols > OMP_THRESHOLD)
            {
                Parallel.For(1, rows - 1, r =>
                {
                    int offset = r * cols;
                    for (int c = 1; c < cols - 1; c++)
                    {
                        int idx = offset + c;
                        res[idx] = x[idx - cols] + x[idx + cols] + x[idx - 1] + x[idx + 1] - 4f * x[idx];
                    }
                });
            }
            else
            {
                for (int r = 1; r < rows - 1; r++)
                {
                    int offset = r * cols;
                    for (int c = 1; c < cols - 1; c++)
                    {
                        int idx = offset + c;
                        res[idx] = x[idx - cols] + x[idx + cols] + x[idx - 1] + x[idx + 1] - 4f * x[idx];
                    }
                }
            }

            // 좌우 경계
            for (int r = 1; r < rows - 1; r++)
            {
                int rOffset = r * cols;
                res[rOffset] = x[rOffset - cols] + x[rOffset + cols] + x[rOffset + 1] - 3f * x[rOffset];
                
                int rEndIdx = rOffset + cols - 1;
                res[rEndIdx] = x[rEndIdx - cols] + x[rEndIdx + cols] + x[rEndIdx - 1] - 3f * x[rEndIdx];
            }

            // 상하 경계
            for (int c = 1; c < cols - 1; c++)
            {
                res[c] = x[cols + c] + x[c - 1] + x[c + 1] - 3f * x[c];
                
                int bIdx = (rows - 1) * cols + c;
                res[bIdx] = x[bIdx - cols] + x[bIdx - 1] + x[bIdx + 1] - 3f * x[bIdx];
            }

            // 4개 코너 점
            res[0] = x[cols] + x[1] - 2f * x[0];
            
            int blIdx = (rows - 1) * cols;
            res[blIdx] = x[blIdx - cols] + x[blIdx + 1] - 2f * x[blIdx];
            
            int trIdx = cols - 1;
            res[trIdx] = x[trIdx + cols] + x[trIdx - 1] - 2f * x[trIdx];
            
            int brIdx = (rows - 1) * cols + cols - 1;
            res[brIdx] = x[brIdx - cols] + x[brIdx - 1] - 2f * x[brIdx];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float Snrm(int n, float[] sx)
        {
            return MathF.Sqrt(DotProduct(sx, sx, n));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void VAdds(float[] a, float s, float[] b, float[] c, int n)
        {
            int i = 0;
            int vectorLength = Vector<float>.Count;
            if (n >= vectorLength)
            {
                var vS = new Vector<float>(s);
                for (; i <= n - vectorLength; i += vectorLength)
                {
                    (new Vector<float>(a, i) + vS * new Vector<float>(b, i)).CopyTo(c, i);
                }
            }
            for (; i < n; i++) c[i] = a[i] + s * b[i];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void VSubs(float[] a, float s, float[] b, float[] c, int n)
        {
            int i = 0;
            int vectorLength = Vector<float>.Count;
            if (n >= vectorLength)
            {
                var vS = new Vector<float>(s);
                for (; i <= n - vectorLength; i += vectorLength)
                {
                    (new Vector<float>(a, i) - vS * new Vector<float>(b, i)).CopyTo(c, i);
                }
            }
            for (; i < n; i++) c[i] = a[i] - s * b[i];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float DotProduct(float[] a, float[] b, int n)
        {
            float sum = 0f;
            int i = 0;
            int vectorLength = Vector<float>.Count;
            if (n >= vectorLength)
            {
                for (; i <= n - vectorLength; i += vectorLength)
                {
                    sum += Vector.Dot(new Vector<float>(a, i), new Vector<float>(b, i));
                }
            }
            for (; i < n; i++) sum += a[i] * b[i];
            return sum;
        }

        /// <summary>
        /// Biconjugate Gradient Method 코어 엔진입니다.
        /// GC 부하를 막기 위해 내부 배열을 ArrayPool에서 빌려 최적화했습니다.
        /// </summary>
        private static void LinBcg(int n, float[] b, float[] x, float tol, int itmax, out int iter, out float err, int rows, int cols)
        {
            var pool = ArrayPool<float>.Shared;
            float[] p = pool.Rent(n);
            float[] pp = pool.Rent(n);
            float[] r = pool.Rent(n);
            float[] rr = pool.Rent(n);
            float[] z = pool.Rent(n);
            float[] zz = pool.Rent(n);

            try
            {
                iter = 0;
                ATimes(x, r, rows, cols);
                
                VSubs(b, 1.0f, r, r, n); // r = b - r
                Array.Copy(r, rr, n);    // rr = r

                ATimes(r, rr, rows, cols); 
                float bnrm = Snrm(n, b);
                ASolve(r, z, n);

                float bkden = 1.0f;
                err = 0f;

                while (iter <= itmax)
                {
                    ++iter;
                    ASolve(rr, zz, n);
                    
                    float bknum = DotProduct(z, rr, n);

                    if (iter == 1)
                    {
                        Array.Copy(z, p, n);
                        Array.Copy(zz, pp, n);
                    }
                    else
                    {
                        float bk = bknum / bkden;
                        VAdds(z, bk, p, p, n);
                        VAdds(zz, bk, pp, pp, n);
                    }
                    
                    bkden = bknum;
                    ATimes(p, z, rows, cols);
                    
                    float akden = DotProduct(z, pp, n);
                    float ak = bknum / akden;
                    
                    ATimes(pp, zz, rows, cols);
                    
                    VAdds(x, ak, p, x, n);
                    VSubs(r, ak, z, r, n);
                    VSubs(rr, ak, zz, rr, n);
                    
                    ASolve(r, z, n);
                    
                    err = Snrm(n, r) / bnrm;
                    if (err <= tol) break;
                }
            }
            finally
            {
                pool.Return(p);
                pool.Return(pp);
                pool.Return(r);
                pool.Return(rr);
                pool.Return(z);
                pool.Return(zz);
            }
        }
    }
}