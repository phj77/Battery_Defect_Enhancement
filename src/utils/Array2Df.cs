using System;
using System.Runtime.CompilerServices;

namespace PdeSolver.Common
{
    /// <summary>
    /// 1차원 배열을 사용하여 2차원 데이터를 표현하는 고성능 배열 클래스.
    /// Multigrid 솔버와 FFT 솔버 모두에서 호환되도록 프로퍼티와 메서드를 통합했습니다.
    /// </summary>
    public class Array2Df
    {
        // Multigrid 솔버 호환용 프로퍼티
        public int Cols { get; }
        public int Rows { get; }

        // FFT 솔버 호환용 프로퍼티 (Cols, Rows를 반환하도록 연결)
        public int Width => Cols;
        public int Height => Rows;

        public float[] Data { get; }
        public int Length => Data.Length;

        public Array2Df(int widthOrCols, int heightOrRows)
        {
            Cols = widthOrCols;
            Rows = heightOrRows;
            Data = new float[Cols * Rows];
        }

        public ref float this[int x, int y]
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => ref Data[y * Cols + x];
        }

        // FFT 솔버 호환용 메서드
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public float Get(int i) => Data[i];

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Set(int i, float value) => Data[i] = value;

        // Multigrid 솔버 호환용 메서드
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Reset() => Array.Clear(Data);
    }
}