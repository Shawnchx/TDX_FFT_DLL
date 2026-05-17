// FFTAnalysis.h
// 傅里叶分析模块：FFT计算、Top9周期提取、正弦拟合
// 通达信DLL插件 - 傅里叶周期分析指标
//
// 优化说明（v5）:
//   1. FFT复数缓冲区改为静态全局，消除每次堆分配/释放
//   2. 旋转因子表预计算，FFT蝶形运算不再调用 cos/sin
//   3. 分量重建用余弦基，增量旋转法避免逐点 cos 调用
//      公式: component(t) = A * cos(ω*t + φ)
//   4. 缓存校验改为 FNV-1a 哈希（O(N)单次，避免重复比较）
//   5. 修复 MapNormSine 对 DataLen 而非 N_USE 遍历的潜在越界
//   6. TOP_K 边界保护，ampBins 不足时安全填充零分量
//   7. (v5) 相位定义改为 t=0 在最新bar（右端）:
//      φ_right = φ_fft + ω*(N_USE-1)，归一化到 [-π, π]
//      物理意义: φ 直接表示"当前bar在该周期的余弦相位角"
//      可横向比较不同股票的相位超前/滞后关系
//      注: 曲线重建代码不变（左端起点初始值等价），仅相位输出值改变

#pragma once
#include <complex>
#include <cmath>
#include <cstring>
#include <algorithm>

// ============================================================
//  常量定义
// ============================================================
static const int FFT_N  = 8192;   // FFT最大点数上限（2^13），实际N_USE动态向下取2次幂
static const int TOP_K  = 9;      // 保留的主要周期数量（Top1~Top9）
static const double PI2 = 6.283185307179586;

// ============================================================
//  存储一个频率分量的全部信息
// ============================================================
struct FreqComponent {
    double frequency;
    double period;
    double amplitude;
    double phase;
    int    binIndex;
};

// ============================================================
//  FFT旋转因子表（全局，程序启动时一次性初始化）
//  twiddle[k] = exp(-j * 2π * k / FFT_N)，只存正变换
//  FFT_N/2 个复数 = 64KB，永久占用但消除所有蝶形运算中的三角计算
// ============================================================
struct TwiddleTable {
    double re[FFT_N / 2];
    double im[FFT_N / 2];
    bool   ready;

    TwiddleTable() : ready(false) {}

    void init() {
        if (ready) return;
        for (int k = 0; k < FFT_N / 2; k++) {
            double ang = -PI2 * k / FFT_N;
            re[k] = cos(ang);
            im[k] = sin(ang);
        }
        ready = true;
    }
};
extern TwiddleTable g_twiddle;

// ============================================================
//  全局缓存
// ============================================================
struct FFTCache {
    int           dataLen;
    int           N_USE;
    uint32_t      inputHash;                    // FNV-1a哈希，用于快速缓存失效检测
    FreqComponent top[TOP_K];
    float         sineCurves[TOP_K][FFT_N];
    float         fitCurve[FFT_N];
    float         trendCurve[FFT_N];
    bool          valid;

    FFTCache() : dataLen(0), N_USE(0), inputHash(0), valid(false) {
        memset(sineCurves,  0, sizeof(sineCurves));
        memset(fitCurve,    0, sizeof(fitCurve));
        memset(trendCurve,  0, sizeof(trendCurve));
        memset(top,         0, sizeof(top));
    }
};
extern FFTCache g_fftCache;

// ============================================================
//  静态复数工作缓冲区（避免每次FFT堆分配）
//  FFT_N个complex<double> = 8192 * 16 = 128KB，全局静态
// ============================================================
extern std::complex<double> g_fftBuf[FFT_N];

// ============================================================
//  FNV-1a 哈希（快速，适合float数组指纹）
// ============================================================
inline uint32_t fnv1a_floats(const float* data, int n)
{
    uint32_t h = 2166136261u;
    const uint8_t* p = reinterpret_cast<const uint8_t*>(data);
    int bytes = n * (int)sizeof(float);
    for (int i = 0; i < bytes; i++) {
        h ^= p[i];
        h *= 16777619u;
    }
    return h;
}

// ============================================================
//  就地 Cooley-Tukey FFT（使用旋转因子表，N必须为2的幂）
// ============================================================
inline void fft_inplace(std::complex<double>* a, int n)
{
    // Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double tr = a[i].real(), ti = a[i].imag();
            a[i] = a[j];
            a[j] = std::complex<double>(tr, ti);
        }
    }
    // Butterfly（使用旋转因子表）
    for (int len = 2; len <= n; len <<= 1) {
        int step = FFT_N / len;  // 旋转因子步长（基于FFT_N归一化的表）
        for (int i = 0; i < n; i += len) {
            for (int j = 0; j < len / 2; j++) {
                int   ki = j * step;          // 查表索引
                double wr = g_twiddle.re[ki];
                double wi = g_twiddle.im[ki];
                double ur = a[i+j].real(),       ui = a[i+j].imag();
                double vr = a[i+j+len/2].real(), vi = a[i+j+len/2].imag();
                double xr = wr*vr - wi*vi,       xi = wr*vi + wi*vr;
                a[i+j]         = std::complex<double>(ur+xr, ui+xi);
                a[i+j+len/2]   = std::complex<double>(ur-xr, ui-xi);
            }
        }
    }
}

// ============================================================
//  主分析函数
// ============================================================
inline void RunFFTAnalysis(int dataLen, const float* pData)
{
    g_twiddle.init();  // 旋转因子表（第一次调用时初始化，后续直接返回）

    // --- 确定 N_USE ---
    int useLen = (dataLen >= FFT_N) ? FFT_N : dataLen;
    int N_USE = 64;
    while (N_USE * 2 <= useLen) N_USE *= 2;

    int srcStart = dataLen - N_USE;
    const float* src = pData + srcStart;

    // --- 哈希校验缓存（O(N)但只做一次，不重复比较）---
    uint32_t h = fnv1a_floats(src, N_USE);
    if (g_fftCache.valid &&
        g_fftCache.dataLen == dataLen &&
        g_fftCache.N_USE   == N_USE   &&
        g_fftCache.inputHash == h)
    {
        return;  // 缓存命中，直接返回
    }

    // --- 更新缓存标识 ---
    g_fftCache.dataLen   = dataLen;
    g_fftCache.N_USE     = N_USE;
    g_fftCache.inputHash = h;
    g_fftCache.valid     = false;  // 计算期间标记无效

    // --- 去线性趋势（最小二乘 y = a*t + b）---
    double sumT=0, sumV=0, sumTT=0, sumTV=0;
    double dN = (double)N_USE;
    for (int i = 0; i < N_USE; i++) {
        double v = src[i];
        sumT  += i;
        sumV  += v;
        sumTT += (double)i * i;
        sumTV += (double)i * v;
    }
    double denom  = dN * sumTT - sumT * sumT;
    double trend_a = (denom != 0.0) ? (dN * sumTV - sumT * sumV) / denom : 0.0;
    double trend_b = (sumV - trend_a * sumT) / dN;

    // --- 填充静态缓冲区（去趋势）---
    for (int i = 0; i < N_USE; i++)
        g_fftBuf[i] = std::complex<double>(src[i] - (trend_a * i + trend_b), 0.0);

    // --- FFT ---
    fft_inplace(g_fftBuf, N_USE);

    // --- 提取幅值并找 Top K ---
    int halfN = N_USE / 2;
    // 使用固定大小的局部数组（最大 FFT_N/2 = 4096），避免vector堆分配
    static double  s_amps[FFT_N / 2];
    static int     s_bins[FFT_N / 2];
    int nBins = halfN - 1;  // bin 1 ~ halfN-1

    for (int k = 1; k < halfN; k++) {
        double re = g_fftBuf[k].real(), im = g_fftBuf[k].imag();
        s_amps[k-1] = sqrt(re*re + im*im) * 2.0 / N_USE;
        s_bins[k-1] = k;
    }

    // 部分排序：只取 Top K（用 partial_sort 替代全排序）
    // 创建索引数组，按幅值降序取前 TOP_K 个
    static int s_idx[FFT_N / 2];
    int actualK = (nBins >= TOP_K) ? TOP_K : nBins;
    for (int i = 0; i < nBins; i++) s_idx[i] = i;
    std::partial_sort(s_idx, s_idx + actualK, s_idx + nBins,
        [](int a, int b){ return s_amps[a] > s_amps[b]; });

    // --- 填充 Top K 分量 ---
    //
    // v5: 相位转换为 t=0 在最新bar（右端）的余弦相位
    //   φ_fft   = atan2(im, re)，对应 t=0 在左端（有效段第一根bar）
    //   φ_right = φ_fft + ω*(N_USE-1)，对应 t=0 在最新bar
    //   物理意义: φ_right 直接表示"当前bar（最新bar）在该周期的余弦相位角"
    //             cos(φ_right) = 当前bar的归一化分量值
    //   归一化到 [-π, π] 方便解读（正值=余弦前半段，负值=后半段）
    //
    for (int r = 0; r < TOP_K; r++) {
        if (r < actualK) {
            int    idx   = s_idx[r];
            int    k     = s_bins[idx];
            double amp   = s_amps[idx];
            double re    = g_fftBuf[k].real(), im = g_fftBuf[k].imag();
            double omega = PI2 * k / N_USE;
            double phi_fft   = atan2(im, re);
            // 转换：最新bar = 旧坐标 t = N_USE-1
            double phi_right = phi_fft + omega * (N_USE - 1);
            // 归一化到 [-π, π]
            while (phi_right >  PI2 / 2) phi_right -= PI2;
            while (phi_right < -PI2 / 2) phi_right += PI2;
            g_fftCache.top[r].binIndex  = k;
            g_fftCache.top[r].frequency = (double)k / N_USE;
            g_fftCache.top[r].period    = (double)N_USE / k;
            g_fftCache.top[r].amplitude = amp;
            g_fftCache.top[r].phase     = phi_right;
        } else {
            // 数据不足时安全填充零分量（避免越界输出垃圾）
            g_fftCache.top[r] = FreqComponent{0.0, 0.0, 0.0, 0.0, 0};
        }
    }

    // --- 生成频域分量曲线（增量旋转法，避免逐点 cos/sin 调用）---
    //
    //  重建公式（余弦基，t=0 定义在最新bar右端）:
    //    component(t) = A * cos(ω*t + φ_right)
    //    其中 t ∈ [-(N_USE-1), 0]，最新bar对应 t=0
    //
    //  存储时从左端开始写（t = -(N_USE-1) 到 t = 0）：
    //    左端 t = -(N_USE-1) 时：
    //      cv_init = cos(ω*(-(N_USE-1)) + φ_right)
    //              = cos(-ω*(N_USE-1) + φ_fft + ω*(N_USE-1))
    //              = cos(φ_fft)          ← 与 v4 完全相同！
    //    所以曲线重建的初始状态 (cv,sv) 不变，仅相位输出含义改变。
    //
    //  增量旋转每步向右（时间增大方向）旋转 +ω。
    //
    int pad = FFT_N - N_USE;

    for (int r = 0; r < TOP_K; r++) {
        // 前 pad 个位置清零（有效数据之前填0）
        memset(g_fftCache.sineCurves[r], 0, pad * sizeof(float));

        if (g_fftCache.top[r].amplitude < 1e-15) {
            // 零幅值分量直接全清零
            memset(g_fftCache.sineCurves[r] + pad, 0, N_USE * sizeof(float));
            continue;
        }

        int    k     = g_fftCache.top[r].binIndex;
        double amp   = g_fftCache.top[r].amplitude;
        double phi   = g_fftCache.top[r].phase;   // φ_right（t=0在最新bar）
        double omega = PI2 * k / N_USE;

        // 增量旋转子（每步旋转 +omega 弧度，时间向右）
        double cosW = cos(omega), sinW = sin(omega);
        // 初始状态对应左端 t=-(N_USE-1):
        //   cv_init = cos(φ_right - ω*(N_USE-1)) = cos(φ_fft)
        double phi_left = phi - omega * (N_USE - 1);
        double cv = cos(phi_left), sv = sin(phi_left);

        for (int t = 0; t < N_USE; t++) {
            g_fftCache.sineCurves[r][pad + t] = (float)(amp * cv);
            // 增量旋转: cos(ω(t+1)+φ) = cv*cosW - sv*sinW
            //           sin(ω(t+1)+φ) = sv*cosW + cv*sinW
            double cv2 = cv * cosW - sv * sinW;
            double sv2 = sv * cosW + cv * sinW;
            cv = cv2; sv = sv2;
        }
    }

    // --- 生成趋势线和全合成曲线 ---
    memset(g_fftCache.fitCurve,   0, pad * sizeof(float));
    memset(g_fftCache.trendCurve, 0, pad * sizeof(float));

    for (int t = 0; t < N_USE; t++) {
        float trend = (float)(trend_a * t + trend_b);
        g_fftCache.trendCurve[pad + t] = trend;
        float val = trend;
        for (int r = 0; r < TOP_K; r++)
            val += g_fftCache.sineCurves[r][pad + t];
        g_fftCache.fitCurve[pad + t] = val;
    }

    g_fftCache.valid = true;
}

// ============================================================
//  辅助：将 FFTCache 序列映射到通达信输出数组
// ============================================================
inline void MapFFTSeqToOutput(int dataLen, float* pfOUT, const float* fftSeq)
{
    int N_USE  = g_fftCache.N_USE;
    int pad    = FFT_N - N_USE;
    int kStart = dataLen - N_USE;

    if (kStart > 0)
        memset(pfOUT, 0, kStart * sizeof(float));

    int copyLen = dataLen - (kStart > 0 ? kStart : 0);
    memcpy(pfOUT + (kStart > 0 ? kStart : 0),
           fftSeq + pad + (kStart > 0 ? 0 : -kStart),
           copyLen * sizeof(float));
}
