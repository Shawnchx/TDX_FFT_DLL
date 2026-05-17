"""
验证 FFT 重建逻辑是否正确
模拟 DLL 的 RunFFTAnalysis 全流程，打印每步的误差
"""
import numpy as np

N = 1024
TOP_K = 6

# ---- 生成一段模拟价格序列（正弦叠加 + 噪声，模拟真实股价）----
np.random.seed(42)
t = np.arange(N)
# 用几个已知周期叠加，理论上FFT应该能完美还原
price = (100
         + 5 * np.cos(2*np.pi*4/N * t + 0.3)   # 周期256
         + 3 * np.cos(2*np.pi*10/N * t + 1.2)  # 周期102
         + 2 * np.cos(2*np.pi*32/N * t - 0.5)  # 周期32
         + 1.5 * np.cos(2*np.pi*64/N * t + 0.8)# 周期16
         + 1 * np.cos(2*np.pi*128/N * t - 1.0) # 周期8
         + 0.8 * np.cos(2*np.pi*8/N * t + 0.5) # 周期128
         + np.random.randn(N) * 0.5)            # 噪声

print(f"输入序列长度: {N}, 均值: {price.mean():.4f}")

# ============================================================
# 步骤1: 去均值
# ============================================================
mean = price.mean()
x = price - mean
print(f"\n[步骤1] 去均值，mean={mean:.4f}")

# ============================================================
# 步骤2: FFT
# ============================================================
X = np.fft.fft(x)
print(f"[步骤2] FFT完成，bin数={len(X)}")

# ============================================================
# 步骤3: 计算幅值（双边归一化 *2/N），取Top6
# ============================================================
half = N // 2
amps = np.abs(X[1:half]) * 2.0 / N
bins = np.argsort(amps)[::-1][:TOP_K] + 1  # +1因为从bin1开始

print(f"\n[步骤3] Top{TOP_K} bin索引: {bins}")
print(f"        对应周期(bar): {[f'{N/k:.1f}' for k in bins]}")
print(f"        对应幅值: {[f'{amps[k-1]:.4f}' for k in bins]}")

# ============================================================
# 步骤4A: DLL当前做法 —— 用sin重建（相位+π/2）
# ============================================================
fit_sin = np.zeros(N) + mean
for k in bins:
    amp = np.abs(X[k]) * 2.0 / N
    phi = np.angle(X[k])          # FFT相位（余弦基）
    phi_sin = phi + np.pi/2       # +π/2 转为sin基
    fit_sin += amp * np.sin(2*np.pi*k/N * t + phi_sin)

rmse_sin = np.sqrt(np.mean((fit_sin - price)**2))
print(f"\n[步骤4A] sin重建（+π/2相位）RMSE = {rmse_sin:.6f}")

# ============================================================
# 步骤4B: 标准做法 —— 直接用cos重建（无相位偏移）
# ============================================================
fit_cos = np.zeros(N) + mean
for k in bins:
    amp = np.abs(X[k]) * 2.0 / N
    phi = np.angle(X[k])
    fit_cos += amp * np.cos(2*np.pi*k/N * t + phi)

rmse_cos = np.sqrt(np.mean((fit_cos - price)**2))
print(f"[步骤4B] cos重建（标准FFT逆变换）RMSE = {rmse_cos:.6f}")

# ============================================================
# 步骤4C: 完整IFFT（所有频率）
# ============================================================
fit_full = np.fft.ifft(X).real + mean
rmse_full = np.sqrt(np.mean((fit_full - price)**2))
print(f"[步骤4C] 完整IFFT（512个频率）RMSE = {rmse_full:.10f}（理论应≈0）")

# ============================================================
# 步骤4D: Top6复数IFFT（最精确的Top6重建）
# ============================================================
X_top6 = np.zeros(N, dtype=complex)
for k in bins:
    X_top6[k] = X[k]
    X_top6[N-k] = X[N-k]  # 共轭对称
fit_top6_ifft = np.fft.ifft(X_top6).real + mean
rmse_top6_ifft = np.sqrt(np.mean((fit_top6_ifft - price)**2))
print(f"[步骤4D] Top6 IFFT重建 RMSE = {rmse_top6_ifft:.6f}（与4A/4B对比）")

# ============================================================
# 诊断：sin vs cos 差异来源
# ============================================================
print(f"\n=== 诊断 ===")
print(f"4A(sin+π/2) vs 4B(cos) 最大差: {np.max(np.abs(fit_sin - fit_cos)):.2e}")
print(f"4B(cos) vs 4D(Top6 IFFT) 最大差: {np.max(np.abs(fit_cos - fit_top6_ifft)):.2e}")

# ============================================================
# 步骤5: 检查补零影响
# ============================================================
print(f"\n=== 补零影响分析 ===")
# 只有500根K线时，前端补0到1024
price_short = price[524:]  # 取后500根
input_padded = np.zeros(N)
input_padded[N-500:] = price_short

mean_padded = input_padded[input_padded != 0].mean()
x_padded = np.where(input_padded != 0, input_padded - mean_padded, 0.0)
X_padded = np.fft.fft(x_padded)
amps_padded = np.abs(X_padded[1:half]) * 2.0 / N
bins_padded = np.argsort(amps_padded)[::-1][:TOP_K] + 1

print(f"补零后 Top{TOP_K} bin: {bins_padded}")
print(f"原始序列 Top{TOP_K} bin: {bins}")
print(f"补零导致的频率识别变化: {'有差异!' if list(bins_padded) != list(bins) else '一致'}")

# 补零后重建误差（仅对500根有数据的部分评估）
fit_padded = np.zeros(N) + mean_padded
t_full = np.arange(N)
for k in bins_padded:
    amp = np.abs(X_padded[k]) * 2.0 / N
    phi = np.angle(X_padded[k])
    fit_padded += amp * np.cos(2*np.pi*k/N * t_full + phi)

rmse_padded = np.sqrt(np.mean((fit_padded[N-500:] - price_short)**2))
print(f"补零序列重建误差(后500根): {rmse_padded:.4f}")
print(f"完整1024序列重建误差: {rmse_cos:.4f}")
print(f"结论: 补零{'严重' if rmse_padded > rmse_cos * 3 else '轻微'}影响拟合精度")
