// TDX_FFT_DLL.cpp
// 通达信傅里叶周期分析插件 DLL
//
// 接口规范：
//   通达信不通过函数名调用，而是调用 RegisterTdxFunc 获取函数注册表。
//   必须导出: BOOL RegisterTdxFunc(PluginTCalcFuncInfo** pFun)
//   部署目录: <通达信安装目录>\T0002\dlls\
//   公式调用: TDXDLL{绑定序号}(内部编号, 输入A, 输入B, 输入C)
//
// ============================================================
//  函数编号规划
// ============================================================
//  1~10  主计算 & 合成拟合曲线
//    1   主计算（触发FFT，直通Close输出）
//    2   纯趋势线（线性拟合 a*t+b，绝对价格）
//    3   趋势 + Top1          合成
//    4   趋势 + Top1~2        合成
//    5   趋势 + Top1~3        合成
//    6   趋势 + Top1~4        合成
//    7   趋势 + Top1~5        合成
//    8   趋势 + Top1~6        合成
//    9   趋势 + Top1~7        合成
//   10   趋势 + Top1~9 全合成（最完整，=fitCurve）
//
//  每个 Top 分量占一段 (N0*10+1 ~ N0*10+6)，N0=1~9
//    x1  幅值
//    x2  周期 (bar)
//    x3  相位 (弧度)
//    x4  频率 (cycles/bar)
//    x5  原始正弦曲线（真实幅值，用于主图 BASE+sin 叠加）
//    x6  归一化正弦（幅值=1，用于副图对比相位和周期）
//
//  示例：
//    11  Top1 幅值   12 Top1 周期   13 Top1 相位
//    14  Top1 频率   15 Top1 正弦   16 Top1 归一化正弦
//    21  Top2 幅值   ...            26 Top2 归一化正弦
//    ...
//    91  Top9 幅值   ...            96 Top9 归一化正弦

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include "PluginTCalcFunc.h"
#include "FFTAnalysis.h"

// ============================================================
//  全局变量定义
// ============================================================
FFTCache                g_fftCache;
TwiddleTable            g_twiddle;
std::complex<double>    g_fftBuf[FFT_N];

// ============================================================
//  DLL 入口
// ============================================================
BOOL APIENTRY DllMain(HMODULE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
{
    if (ul_reason_for_call == DLL_PROCESS_ATTACH) {
        g_fftCache = FFTCache();
        g_twiddle.init();   // 预计算旋转因子表，加载时完成，后续FFT无需三角计算
    }
    return TRUE;
}

// ============================================================
//  辅助宏
// ============================================================
// 输出常量序列（所有bar相同值）
#define OUTPUT_CONST(val) \
    for (int i = 0; i < DataLen; i++) pfOUT[i] = (float)(val);

// 确保FFT已计算（以pfINa作为价格输入）
#define ENSURE_FFT() \
    RunFFTAnalysis(DataLen, pfINa); \
    if (!g_fftCache.valid) { OUTPUT_CONST(0.0f); return; }

// ============================================================
//  Func1 - 主计算（触发FFT，直通Close）
// ============================================================
static void Func1(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{
    RunFFTAnalysis(DataLen, pfINa);
    for (int i = 0; i < DataLen; i++) pfOUT[i] = pfINa[i];
}

// ============================================================
//  合成拟合曲线辅助
//  nComp=0: 纯趋势线
//  nComp=1~9: 趋势 + 前nComp个正弦分量
//  算法: full_fit[t] = trend[t] + sum(sine[0..8][t])
//        partial[t]  = full_fit[t] - sum(sine[nComp..8][t])
//                    = trend[t] + sum(sine[0..nComp-1][t])
// ============================================================
static void MapPartialFit(int DataLen, float* pfOUT, int nComp)
{
    int N_USE  = g_fftCache.N_USE;
    int pad    = FFT_N - N_USE;
    int kStart = DataLen - N_USE;

    for (int i = 0; i < DataLen; i++) {
        if (i < kStart) {
            pfOUT[i] = 0.0f;
        } else {
            int t = i - kStart;
            if (nComp == 0) {
                // 纯趋势线
                pfOUT[i] = g_fftCache.trendCurve[pad + t];
            } else {
                // 全合成减去不需要的高阶分量
                float val = g_fftCache.fitCurve[pad + t];
                for (int r = nComp; r < TOP_K; r++)
                    val -= g_fftCache.sineCurves[r][pad + t];
                pfOUT[i] = val;
            }
        }
    }
}

// Func2~10: 合成拟合曲线
static void Func2 (int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapPartialFit(DataLen, pfOUT, 0); }  // 纯趋势线
static void Func3 (int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapPartialFit(DataLen, pfOUT, 1); }  // 趋势+Top1
static void Func4 (int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapPartialFit(DataLen, pfOUT, 2); }  // 趋势+Top1~2
static void Func5 (int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapPartialFit(DataLen, pfOUT, 3); }  // 趋势+Top1~3
static void Func6 (int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapPartialFit(DataLen, pfOUT, 4); }  // 趋势+Top1~4
static void Func7 (int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapPartialFit(DataLen, pfOUT, 5); }  // 趋势+Top1~5
static void Func8 (int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapPartialFit(DataLen, pfOUT, 6); }  // 趋势+Top1~6
static void Func9 (int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapPartialFit(DataLen, pfOUT, 7); }  // 趋势+Top1~7
static void Func10(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.fitCurve); } // 趋势+Top1~9全合成

// ============================================================
//  归一化正弦输出辅助
//  Bug修复：原版对整个 DataLen 除以幅值，前段0值无影响但逻辑不对；
//           修复为先 MapFFTSeqToOutput，再只对有效段 [kStart, DataLen) 归一化
// ============================================================
static void MapNormSine(int DataLen, float* pfOUT, int rank)
{
    float amp = (float)g_fftCache.top[rank].amplitude;
    MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[rank]);
    if (amp < 1e-10f) return;  // 幅值极小时保持0输出，避免除零放大噪声
    int kStart = DataLen - g_fftCache.N_USE;
    if (kStart < 0) kStart = 0;
    float invAmp = 1.0f / amp;
    for (int i = kStart; i < DataLen; i++)
        pfOUT[i] *= invAmp;
}

// ============================================================
//  Top1~Top9 参数和曲线
//  编号规则: rank N0 (0~8) -> 编号 (N0+1)*10+x, x=1~6
//    x1 幅值  x2 周期  x3 相位  x4 频率  x5 原始正弦  x6 归一化正弦
// ============================================================

// --- Top1 (11~16) ---
static void Func11(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[0].amplitude); }
static void Func12(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[0].period);    }
static void Func13(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[0].phase);     }
static void Func14(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[0].frequency); }
static void Func15(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[0]); }
static void Func16(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapNormSine(DataLen, pfOUT, 0); }

// --- Top2 (21~26) ---
static void Func21(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[1].amplitude); }
static void Func22(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[1].period);    }
static void Func23(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[1].phase);     }
static void Func24(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[1].frequency); }
static void Func25(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[1]); }
static void Func26(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapNormSine(DataLen, pfOUT, 1); }

// --- Top3 (31~36) ---
static void Func31(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[2].amplitude); }
static void Func32(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[2].period);    }
static void Func33(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[2].phase);     }
static void Func34(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[2].frequency); }
static void Func35(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[2]); }
static void Func36(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapNormSine(DataLen, pfOUT, 2); }

// --- Top4 (41~46) ---
static void Func41(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[3].amplitude); }
static void Func42(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[3].period);    }
static void Func43(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[3].phase);     }
static void Func44(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[3].frequency); }
static void Func45(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[3]); }
static void Func46(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapNormSine(DataLen, pfOUT, 3); }

// --- Top5 (51~56) ---
static void Func51(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[4].amplitude); }
static void Func52(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[4].period);    }
static void Func53(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[4].phase);     }
static void Func54(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[4].frequency); }
static void Func55(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[4]); }
static void Func56(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapNormSine(DataLen, pfOUT, 4); }

// --- Top6 (61~66) ---
static void Func61(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[5].amplitude); }
static void Func62(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[5].period);    }
static void Func63(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[5].phase);     }
static void Func64(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[5].frequency); }
static void Func65(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[5]); }
static void Func66(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapNormSine(DataLen, pfOUT, 5); }

// --- Top7 (71~76) ---
static void Func71(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[6].amplitude); }
static void Func72(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[6].period);    }
static void Func73(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[6].phase);     }
static void Func74(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[6].frequency); }
static void Func75(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[6]); }
static void Func76(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapNormSine(DataLen, pfOUT, 6); }

// --- Top8 (81~86) ---
static void Func81(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[7].amplitude); }
static void Func82(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[7].period);    }
static void Func83(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[7].phase);     }
static void Func84(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[7].frequency); }
static void Func85(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[7]); }
static void Func86(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapNormSine(DataLen, pfOUT, 7); }

// --- Top9 (91~96) ---
static void Func91(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[8].amplitude); }
static void Func92(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[8].period);    }
static void Func93(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[8].phase);     }
static void Func94(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); OUTPUT_CONST(g_fftCache.top[8].frequency); }
static void Func95(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapFFTSeqToOutput(DataLen, pfOUT, g_fftCache.sineCurves[8]); }
static void Func96(int DataLen, float* pfOUT, float* pfINa, float* pfINb, float* pfINc)
{ ENSURE_FFT(); MapNormSine(DataLen, pfOUT, 8); }

// ============================================================
//  函数注册表
// ============================================================
static PluginTCalcFuncInfo g_CalcFuncSets[] =
{
    // 主计算 & 合成拟合 (1~10)
    {  1, (pPluginFUNC)&Func1  },  // 主计算，直通Close
    {  2, (pPluginFUNC)&Func2  },  // 纯趋势线
    {  3, (pPluginFUNC)&Func3  },  // 趋势+Top1
    {  4, (pPluginFUNC)&Func4  },  // 趋势+Top1~2
    {  5, (pPluginFUNC)&Func5  },  // 趋势+Top1~3
    {  6, (pPluginFUNC)&Func6  },  // 趋势+Top1~4
    {  7, (pPluginFUNC)&Func7  },  // 趋势+Top1~5
    {  8, (pPluginFUNC)&Func8  },  // 趋势+Top1~6
    {  9, (pPluginFUNC)&Func9  },  // 趋势+Top1~7
    { 10, (pPluginFUNC)&Func10 },  // 趋势+Top1~9 全合成
    // Top1 (11~16): 幅值/周期/相位/频率/正弦/归一化
    { 11, (pPluginFUNC)&Func11 }, { 12, (pPluginFUNC)&Func12 },
    { 13, (pPluginFUNC)&Func13 }, { 14, (pPluginFUNC)&Func14 },
    { 15, (pPluginFUNC)&Func15 }, { 16, (pPluginFUNC)&Func16 },
    // Top2 (21~26)
    { 21, (pPluginFUNC)&Func21 }, { 22, (pPluginFUNC)&Func22 },
    { 23, (pPluginFUNC)&Func23 }, { 24, (pPluginFUNC)&Func24 },
    { 25, (pPluginFUNC)&Func25 }, { 26, (pPluginFUNC)&Func26 },
    // Top3 (31~36)
    { 31, (pPluginFUNC)&Func31 }, { 32, (pPluginFUNC)&Func32 },
    { 33, (pPluginFUNC)&Func33 }, { 34, (pPluginFUNC)&Func34 },
    { 35, (pPluginFUNC)&Func35 }, { 36, (pPluginFUNC)&Func36 },
    // Top4 (41~46)
    { 41, (pPluginFUNC)&Func41 }, { 42, (pPluginFUNC)&Func42 },
    { 43, (pPluginFUNC)&Func43 }, { 44, (pPluginFUNC)&Func44 },
    { 45, (pPluginFUNC)&Func45 }, { 46, (pPluginFUNC)&Func46 },
    // Top5 (51~56)
    { 51, (pPluginFUNC)&Func51 }, { 52, (pPluginFUNC)&Func52 },
    { 53, (pPluginFUNC)&Func53 }, { 54, (pPluginFUNC)&Func54 },
    { 55, (pPluginFUNC)&Func55 }, { 56, (pPluginFUNC)&Func56 },
    // Top6 (61~66)
    { 61, (pPluginFUNC)&Func61 }, { 62, (pPluginFUNC)&Func62 },
    { 63, (pPluginFUNC)&Func63 }, { 64, (pPluginFUNC)&Func64 },
    { 65, (pPluginFUNC)&Func65 }, { 66, (pPluginFUNC)&Func66 },
    // Top7 (71~76)
    { 71, (pPluginFUNC)&Func71 }, { 72, (pPluginFUNC)&Func72 },
    { 73, (pPluginFUNC)&Func73 }, { 74, (pPluginFUNC)&Func74 },
    { 75, (pPluginFUNC)&Func75 }, { 76, (pPluginFUNC)&Func76 },
    // Top8 (81~86)
    { 81, (pPluginFUNC)&Func81 }, { 82, (pPluginFUNC)&Func82 },
    { 83, (pPluginFUNC)&Func83 }, { 84, (pPluginFUNC)&Func84 },
    { 85, (pPluginFUNC)&Func85 }, { 86, (pPluginFUNC)&Func86 },
    // Top9 (91~96)
    { 91, (pPluginFUNC)&Func91 }, { 92, (pPluginFUNC)&Func92 },
    { 93, (pPluginFUNC)&Func93 }, { 94, (pPluginFUNC)&Func94 },
    { 95, (pPluginFUNC)&Func95 }, { 96, (pPluginFUNC)&Func96 },
    {  0, NULL }  // 结束标记
};

// ============================================================
//  RegisterTdxFunc - 唯一需要导出的函数
// ============================================================
extern "C" __declspec(dllexport)
BOOL RegisterTdxFunc(PluginTCalcFuncInfo** pFun)
{
    if (*pFun == NULL)
    {
        *pFun = g_CalcFuncSets;
        return TRUE;
    }
    return FALSE;
}
