using System;
using System.Threading.Tasks;
using PdeSolver.Common;

namespace PreProcessing;

/// <summary>
/// Fattal TMO 실행 전, 입력 이미지 데이터에 대한 전처리를 담당합니다.
/// </summary>
public static class PreProcessor
{
    /// <summary>
    /// 입력 이미지의 각 픽셀에 감마 보정을 적용합니다.
    /// 공식: Out = In ^ gamma
    /// </summary>
    /// <param name="image">전처리를 수행할 이미지 데이터 (Array2Df)</param>
    /// <param name="gamma">적용할 감마 값</param>
    public static void ApplyGamma(Array2Df image, float gamma)
    {
        if (MathF.Abs(gamma - 1.0f) < 1e-6f) return; // 감마가 1이면 연산 생략

        int length = image.Length;
        float[] data = image.Data;

        // .NET 7/C# 11의 Parallel.For를 이용한 병렬 처리
        Parallel.For(0, length, i =>
        {
            // 0 이하의 값이 입력될 경우를 대비해 Clamp 처리 후 거듭제곱 연산
            float val = data[i];
            if (val > 0)
            {
                data[i] = MathF.Pow(val, gamma);
            }
            else
            {
                data[i] = 0.0f;
            }
        });
    }

    /// <summary>
    /// 이미지의 수치 범위를 정규화하거나 특정 상수를 곱하는 등의 추가 전처리를 수행할 수 있습니다.
    /// </summary>
    public static void Normalize(Array2Df image, float scale)
    {
        int length = image.Length;
        float[] data = image.Data;

        Parallel.For(0, length, i =>
        {
            data[i] *= scale;
        });
    }
}