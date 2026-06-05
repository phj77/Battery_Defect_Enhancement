using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using FattalToneMapping;
using FattalToneMapping.Config;
using PdeSolver.Common;

using StbImageSharp;
using StbImageWriteSharp;

using System;
using System.Diagnostics; // Stopwatch를 사용하기 위해 필요
using System.Threading;


string basePath = AppDomain.CurrentDomain.BaseDirectory;
string configPath = Path.Combine(basePath, "config.json");

GlobalTimer.Stopwatch.Start();

if (!File.Exists(configPath))
{
    Console.WriteLine($"cannof find config.json. (search path: {configPath})");
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
        // 입력 경로가 폴더인지 파일인지 확인
        string[] targetFiles;
        if (Directory.Exists(sweep.InputPath))
        {
            // 폴더 내 모든 파일을 가져옵니다. 특정 확장자(예: .hdr)만 필요하다면 
            // Directory.GetFiles(sweep.InputPath, "*.hdr") 로 변경하십시오.
            targetFiles = Directory.GetFiles(sweep.InputPath);
        }
        else if (File.Exists(sweep.InputPath))
        {
            // 단일 파일일 경우
            targetFiles = new[] { sweep.InputPath };
        }
        else
        {
            Console.WriteLine($"[warning] path not found: {sweep.InputPath}");
            continue;
        }

        // 대상 파일 각각에 대해 파라미터 조합 생성
        foreach (var file in targetFiles)
        {
            // 파일명 덮어쓰기 방지를 위해 원본 파일명을 추출
            string fileNameWithoutExt = Path.GetFileNameWithoutExtension(file);

            foreach (var a in sweep.Alphas)
            {
                foreach (var b in sweep.Betas)
                {
                    foreach (var s in sweep.Saturations)
                    {
                        // OutputPrefix, 원본 파일명, 파라미터를 조합하여 결과 파일명 생성
                        // 예: prefix_memorial_a0.1_b0.8_s0.6.png
                        string outPrefix = string.IsNullOrEmpty(sweep.OutputPrefix) ? "" : $"{sweep.OutputPrefix}_";
                        string outName = $"{outPrefix}{fileNameWithoutExt}_a{a}_b{b}_s{s}.png";

                        allTasks.Add(new TonemappingTask(
                            InputPath: file, // 순회 중인 개별 파일 경로 사용
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
}

if (allTasks.Count == 0)
{
    Console.WriteLine("no works to execute");
    return;
}

Console.WriteLine($"start {allTasks.Count} number of tonemapping...\n");

// 3. 작업 순차 실행
foreach (var task in allTasks)
{
    try
    {
        Console.WriteLine($"[progressing...] {task.OutputPath} generating... (Alpha:{task.Alpha}, Beta:{task.Beta})");

        // TODO: 실제 HDR 이미지 로드 로직 구현 필요 (StbImageSharp 등 사용)
        if (!LoadHdrImage(task.InputPath, out int width, out int height, out Array2Df r, out Array2Df g, out Array2Df b))
        {
            Console.WriteLine($"  -> fail to load image: {task.InputPath}");
            continue;
        }

        // 알고리즘 실행
        PfsTmoFattal02.Apply(
            r, g, b, 
            task.Alpha, task.Beta, task.Saturation, task.Noise, 
            task.NewFattal, task.FftSolver, task.DetailLevel
        );


        Console.WriteLine($"makingintensity range [0,1] to [0,255] start: {GlobalTimer.ElapsedSeconds:F2}s");
        // 0~1 사이의 Float 값을 0~255 사이의 Byte 픽셀 데이터로 변환
        byte[] pixels = new byte[width * height * 3];
        for (int i = 0; i < width * height; i++)
        {
            pixels[i * 3 + 0] = (byte)Math.Clamp(r.Data[i] * 255f, 0, 255);
            pixels[i * 3 + 1] = (byte)Math.Clamp(g.Data[i] * 255f, 0, 255);
            pixels[i * 3 + 2] = (byte)Math.Clamp(b.Data[i] * 255f, 0, 255);
        }

        // TODO: 결과 이미지(LDR) 저장 로직 구현 필요
        SaveLdrImage(task.OutputPath, width, height, pixels);
        Console.WriteLine($"  -> complete!");
    }
    catch (Exception ex)
    {
        Console.WriteLine($"  -> [error] {ex.Message}");
    }
}

GlobalTimer.Stopwatch.Stop();
double seconds = GlobalTimer.ElapsedSeconds;
Console.WriteLine($"execution time: {seconds:F4} s");
Console.WriteLine("\nevery experiment over.");
Console.WriteLine("\nprogram is over. press any button.");
Console.ReadKey();


// --- 실제 이미지 입출력 함수 ---

bool LoadHdrImage(string path, out int w, out int h, out Array2Df r, out Array2Df g, out Array2Df b)
{
    Console.WriteLine($"image loading start: {GlobalTimer.ElapsedSeconds:F2}s");

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

    Console.WriteLine($"image loaded over at: {GlobalTimer.ElapsedSeconds:F2}s");

    return true;
}

void SaveLdrImage(string path, int w, int h, byte[] pixels)
{
    Console.WriteLine($"image saving start: {GlobalTimer.ElapsedSeconds:F2}s");

    // 출력 폴더가 없다면 생성
    string directory = Path.GetDirectoryName(path);
    if (!string.IsNullOrEmpty(directory) && !Directory.Exists(directory))
    {
        Directory.CreateDirectory(directory);
    }

    // PNG 파일로 저장
    using Stream stream = File.OpenWrite(path);
    ImageWriter writer = new ImageWriter();
    writer.WritePng(pixels, w, h, StbImageWriteSharp.ColorComponents.RedGreenBlue, stream);

    Console.WriteLine($"image saved over at: {GlobalTimer.ElapsedSeconds:F2}s");
}

public static class GlobalTimer
{
    // 프로그램 전역에서 공유할 고정밀 스톱워치
    public static readonly Stopwatch Stopwatch = new Stopwatch();

    // 편리한 사용을 위해 초 단위 누적 시간을 반환하는 프로퍼티
    public static double ElapsedSeconds => Stopwatch.Elapsed.TotalSeconds;
}