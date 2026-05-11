using System;
using System.Threading.Tasks;
using PdeSolver.Common;
using PostProcessing;

namespace FattalToneMapping
{
    /// <summary>
    /// Fattal 2002 모델을 사용하여 RGB 채널의 톤 매핑을 수행하고 후처리를 적용하는 클래스.
    /// </summary>
    public static class PfsTmoFattal02
    {
        private const float Epsilon = 1e-4f;

        /// <summary>
        /// 원본 RGB 채널 데이터에 Fattal 알고리즘을 적용하고 결과를 반환합니다.
        /// 데이터는 원본 배열(r, g, b)에 덮어씌워집니다.
        /// </summary>
        public static void Apply(Array2Df r, Array2Df g, Array2Df b, 
                                 float optAlpha, float optBeta, float optSaturation, float optNoise, 
                                 bool newFattal, bool fftSolver, int detailLevel, 
                                 Action<int> progressCallback = null)
        {
            if (r == null || g == null || b == null)
            {
                throw new ArgumentNullException("RGB 채널 중 하나 이상이 누락되었습니다.");
            }

            int w = r.Cols;
            int h = r.Rows;
            int length = r.Length;

            // FFT 솔버 사용 시 newFattal 파라미터를 강제로 활성화
            if (fftSolver)
            {
                newFattal = true;
            }

            // 노이즈 파라미터가 0 이하일 경우 기본값 설정
            if (optNoise <= 0.0f)
            {
                optNoise = optAlpha * 0.01f;
            }

            progressCallback?.Invoke(0);

            Array2Df Yr = new Array2Df(w, h);
            Array2Df L = new Array2Df(w, h);

            // 1. 휘도(Luminance) 채널 추출 (RGB -> Y 변환)
            // 국제 표준 Rec. 709 휘도 공식을 사용하여 계산
            Parallel.For(0, length, i =>
            {
                Yr.Data[i] = 0.2126f * r.Data[i] + 0.7152f * g.Data[i] + 0.0722f * b.Data[i];
            });

            // 2. 핵심 Fattal 톤 매핑 알고리즘 실행
            try
            {
                TmoFattal02.Process(w, h, Yr, L, optAlpha, optBeta, optNoise, newFattal, fftSolver, detailLevel, progressCallback);
            }
            catch (Exception ex)
            {
                throw new Exception("Tonemapping Failed!", ex);
            }

            // 3. 색상 재구성 (Color Reconstruction)
            // 압축된 휘도(L)와 원본 휘도(Yr)의 비율을 바탕으로 RGB 채널을 복원하고 채도(Saturation)를 조절
            Parallel.For(0, length, i =>
            {
                float yVal = MathF.Max(Yr.Data[i], Epsilon);
                float lVal = MathF.Max(L.Data[i], Epsilon);

                float rRatio = MathF.Max(r.Data[i] / yVal, 0f);
                float gRatio = MathF.Max(g.Data[i] / yVal, 0f);
                float bRatio = MathF.Max(b.Data[i] / yVal, 0f);

                r.Data[i] = MathF.Pow(rRatio, optSaturation) * lVal;
                g.Data[i] = MathF.Pow(gRatio, optSaturation) * lVal;
                b.Data[i] = MathF.Pow(bRatio, optSaturation) * lVal;
            });

            // 4. 후처리 (Post Processing)
            // 톤 매핑된 각 채널에 디스플레이 감마(2.2) 보정을 적용하고, 안전하게 [0, 1] 범위로 클램핑
            PostProcessor.ApplyDisplayGamma(r, 2.2f);
            PostProcessor.ApplyDisplayGamma(g, 2.2f);
            PostProcessor.ApplyDisplayGamma(b, 2.2f);

            PostProcessor.Clamp(r);
            PostProcessor.Clamp(g);
            PostProcessor.Clamp(b);

            progressCallback?.Invoke(100);
        }
    }
}