using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;
using FattalToneMapping;
using FattalToneMapping.Config;
using PdeSolver.Common;

using StbImageSharp;
using StbImageWriteSharp;

string basePath = AppDomain.CurrentDomain.BaseDirectory;
string configPath = Path.Combine(basePath, "config.json");

if (!File.Exists(configPath))
{
    Console.WriteLine($"설정 파일(config.json)을 찾을 수 없습니다. (탐색 경로: {configPath})");
    return;
}

string jsonContent = File.ReadAllText(configPath);
var config = JsonSerializer.Deserialize<AppConfig>(jsonContent, new JsonSerializerOptions { PropertyNameCaseInsensitive = true });

// 실행할 모든 작업을 담을 리스트
var allTasks = new List<TonemappingTask>();

// 1. 단일 작업 추가
if (config?.SingleTasks != null)
{
    allTasks.AddRange(config.SingleTasks);
}

// 2. 실험(Sweep) 작업을 조합하여 단일 작업으로 풀어서 추가
if (config?.SweepTasks != null)
{
    foreach (var sweep in config.SweepTasks)
    {
        foreach (var a in sweep.Alphas)
        {
            foreach (var b in sweep.Betas)
            {
                foreach (var s in sweep.Saturations)
                {
                    // 파라미터 값이 들어간 파일명 자동 생성 (예: memorial_test_a0.1_b0.8_s0.6.png)
                    string outName = $"{sweep.OutputPrefix}_a{a}_b{b}_s{s}.png";
                    
                    allTasks.Add(new TonemappingTask(
                        InputPath: sweep.InputPath,
                        OutputPath: outName,
                        Alpha: a,
                        Beta: b,
                        Saturation: s,
                        Noise: sweep.Noise,
                        NewFattal: sweep.NewFattal,
                        FftSolver: sweep.FftSolver,
                        DetailLevel: sweep.DetailLevel
                    ));
                }
            }
        }
    }
}

if (allTasks.Count == 0)
{
    Console.WriteLine("실행할 작업이 없습니다.");
    return;
}

Console.WriteLine($"총 {allTasks.Count}개의 톤매핑 작업을 시작합니다...\n");

// 3. 작업 순차 실행
foreach (var task in allTasks)
{
    try
    {
        Console.WriteLine($"[진행 중] {task.OutputPath} 생성 중... (Alpha:{task.Alpha}, Beta:{task.Beta})");

        // TODO: 실제 HDR 이미지 로드 로직 구현 필요 (StbImageSharp 등 사용)
        if (!LoadHdrImage(task.InputPath, out int width, out int height, out Array2Df r, out Array2Df g, out Array2Df b))
        {
            Console.WriteLine($"  -> 이미지 로드 실패: {task.InputPath}");
            continue;
        }

        // 알고리즘 실행
        PfsTmoFattal02.Apply(
            r, g, b, 
            task.Alpha, task.Beta, task.Saturation, task.Noise, 
            task.NewFattal, task.FftSolver, task.DetailLevel
        );

        // TODO: 결과 이미지(LDR) 저장 로직 구현 필요
        SaveLdrImage(task.OutputPath, width, height, r, g, b);
        Console.WriteLine($"  -> 완료!");
    }
    catch (Exception ex)
    {
        Console.WriteLine($"  -> [오류 발생] {ex.Message}");
    }
}

Console.WriteLine("\n모든 실험이 종료되었습니다.");


// --- 실제 이미지 입출력 함수 ---

bool LoadHdrImage(string path, out int w, out int h, out Array2Df r, out Array2Df g, out Array2Df b)
{
    w = 0; h = 0; r = null; g = null; b = null;
    
    if (!File.Exists(path)) return false;

    // HDR 이미지를 Float 배열로 읽어들임
    using Stream stream = File.OpenRead(path);
    ImageResultFloat image = ImageResultFloat.FromStream(stream, StbImageSharp.ColorComponents.RedGreenBlue);
    
    w = image.Width;
    h = image.Height;
    r = new Array2Df(w, h);
    g = new Array2Df(w, h);
    b = new Array2Df(w, h);

    // 1차원 배열로 펼쳐진 RGB 데이터를 분리하여 Array2Df에 저장
    for (int i = 0; i < w * h; i++)
    {
        r.Data[i] = image.Data[i * 3 + 0];
        g.Data[i] = image.Data[i * 3 + 1];
        b.Data[i] = image.Data[i * 3 + 2];
    }

    return true;
}

void SaveLdrImage(string path, int w, int h, Array2Df r, Array2Df g, Array2Df b)
{
    // 출력 폴더가 없다면 생성
    string directory = Path.GetDirectoryName(path);
    if (!string.IsNullOrEmpty(directory) && !Directory.Exists(directory))
    {
        Directory.CreateDirectory(directory);
    }

    // 0~1 사이의 Float 값을 0~255 사이의 Byte 픽셀 데이터로 변환
    byte[] pixels = new byte[w * h * 3];
    for (int i = 0; i < w * h; i++)
    {
        pixels[i * 3 + 0] = (byte)Math.Clamp(r.Data[i] * 255f, 0, 255);
        pixels[i * 3 + 1] = (byte)Math.Clamp(g.Data[i] * 255f, 0, 255);
        pixels[i * 3 + 2] = (byte)Math.Clamp(b.Data[i] * 255f, 0, 255);
    }

    // PNG 파일로 저장
    using Stream stream = File.OpenWrite(path);
    ImageWriter writer = new ImageWriter();
    writer.WritePng(pixels, w, h, StbImageWriteSharp.ColorComponents.RedGreenBlue, stream);
}