using System;
using System.Threading.Tasks;
using PdeSolver.Common;

namespace PostProcessing;

/// <summary>
/// Fattal TMO 실행 후, 결과 이미지 데이터에 대한 후처리를 담당합니다.
/// </summary>
public static class PostProcessor
{
    /// <summary>
    /// TMO 결과물에 최종 감마 보정을 적용하여 디스플레이에 적합한 상태로 만듭니다.
    /// 공식: Out = In ^ (1/gamma)
    /// </summary>
    /// <param name="image">후처리를 수행할 이미지 데이터 (Array2Df)</param>
    /// <param name="gamma">디스플레이 감마 값 (보통 2.2)</param>
    public static void ApplyDisplayGamma(Array2Df image, float gamma)
    {
        if (MathF.Abs(gamma - 1.0f) < 1e-6f) return;

        float invGamma = 1.0f / gamma;
        int length = image.Length;
        float[] data = image.Data;

        Parallel.For(0, length, i =>
        {
            float val = data[i];
            // TMO 결과는 0~1 사이로 클램핑되어 있으나 안전을 위해 다시 확인
            val = Math.Clamp(val, 0.0f, 1.0f);
            data[i] = MathF.Pow(val, invGamma);
        });
    }

    /// <summary>
    /// 최종 출력 전 픽셀 값이 [0, 1] 범위를 벗어나지 않도록 강제합니다.
    /// </summary>
    public static void Clamp(Array2Df image)
    {
        int length = image.Length;
        float[] data = image.Data;

        Parallel.For(0, length, i =>
        {
            data[i] = Math.Clamp(data[i], 0.0f, 1.0f);
        });
    }
}