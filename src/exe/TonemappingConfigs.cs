using System.Collections.Generic;

namespace FattalToneMapping.Config;

// 단일 작업용 (기존과 동일)
public record TonemappingTask(
    string InputPath,
    string OutputPath,
    float Alpha = 0.1f,
    float Beta = 0.8f,
    float Saturation = 0.5f,
    float Noise = 0.001f,
    bool NewFattal = true,
    bool FftSolver = false,
    int DetailLevel = 0
);

// 범위 실험(Parameter Sweep)용 작업
public class SweepTask
{
    public string InputPath { get; set; } = "";
    public string OutputPrefix { get; set; } = "result"; // 파일명 앞에 붙을 이름
    
    // 여러 개의 값을 배열로 넣으면, 모든 조합을 알아서 테스트함
    public List<float> Alphas { get; set; } = new() { 0.1f };
    public List<float> Betas { get; set; } = new() { 0.8f };
    public List<float> Saturations { get; set; } = new() { 0.5f };
    
    // 이 값들은 고정
    public float Noise { get; set; } = 0.001f;
    public bool NewFattal { get; set; } = true;
    public bool FftSolver { get; set; } = false;
    public int DetailLevel { get; set; } = 0;
}

// 전체 설정
public class AppConfig
{
    public List<TonemappingTask> SingleTasks { get; set; } = new();
    public List<SweepTask> SweepTasks { get; set; } = new();
}