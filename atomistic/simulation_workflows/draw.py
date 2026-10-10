import numpy as np
import matplotlib.pyplot as plt

# ========== 参数设置，可按需要调整 ==========
data_file = "stress_strain.dat"

# 平滑窗口（奇数），点数越多越平滑
smooth_window = 21

# 认为是“锯齿事件”的应力跌落阈值（GPa），建议先试 1.0
serration_drop_threshold = 1.0

# 拟合弹性模量时使用的最大工程应变
elastic_fit_max_strain = 0.01   # 1%
# =========================================


def moving_average(y, window):
    """
    简单滑动平均，用于平滑曲线。
    window 需要是奇数。
    """
    if window < 2:
        return y
    window = int(window)
    if window % 2 == 0:
        window += 1
    kernel = np.ones(window) / window
    y_padded = np.pad(y, (window // 2, window // 2), mode="edge")
    y_smooth = np.convolve(y_padded, kernel, mode="valid")
    return y_smooth


# 1. 读取数据
data = np.loadtxt(data_file)
strain_eng = data[:, 0]
stress_eng = data[:, 3]   # sigma_zz (GPa)

# 2. 平滑工程应力
stress_eng_smooth = moving_average(stress_eng, smooth_window)

# 3. 计算真应力–真应变（假设体积不变）
strain_true = np.log(1.0 + strain_eng)
stress_true = stress_eng * (1.0 + strain_eng)
stress_true_smooth = moving_average(stress_true, smooth_window)

# 4. 自动计算初始弹性模量（对工程应力–应变做线性拟合）
mask_elastic = strain_eng <= elastic_fit_max_strain
coeffs = np.polyfit(strain_eng[mask_elastic], stress_eng[mask_elastic], 1)
E_fit = coeffs[0]   # GPa
print(f"Estimated elastic modulus (fit up to {elastic_fit_max_strain:.3f} strain): {E_fit:.2f} GPa")

# 5. 屈服点：平滑曲线上的最大值
yield_index = np.argmax(stress_eng_smooth)
yield_strain = strain_eng[yield_index]
yield_stress = stress_eng_smooth[yield_index]
print(f"Yield point (by max smoothed stress): strain = {yield_strain:.4f}, stress = {yield_stress:.2f} GPa")

# 6. 锯齿事件：从局部峰值到紧接着的谷值的“掉落”超过阈值就算一次 serration
drop_min = serration_drop_threshold  # 单位 GPa

s = stress_eng_smooth
n = len(s)

# 找局部最大值（峰）和最小值（谷）
local_max = []
local_min = []
for i in range(1, n-1):
    if s[i] > s[i-1] and s[i] >= s[i+1]:
        local_max.append(i)
    if s[i] < s[i-1] and s[i] <= s[i+1]:
        local_min.append(i)
local_max = np.array(local_max, dtype=int)
local_min = np.array(local_min, dtype=int)

serration_max_idx = []
serration_min_idx = []

# 对每个局部峰，找它后面的第一个局部谷，看 drop 是否超过阈值
for max_i in local_max:
    mins_after = local_min[local_min > max_i]
    if len(mins_after) == 0:
        continue
    min_i = mins_after[0]
    drop = s[max_i] - s[min_i]
    if drop >= drop_min:
        serration_max_idx.append(max_i)
        serration_min_idx.append(min_i)

serration_max_idx = np.array(serration_max_idx, dtype=int)
serration_min_idx = np.array(serration_min_idx, dtype=int)

print(f"Detected serration events (peak-to-valley drop >= {drop_min} GPa): {len(serration_max_idx)} events")
for k, (imax, imin) in enumerate(zip(serration_max_idx, serration_min_idx), 1):
    print(f"  Event {k}: peak at strain={strain_eng[imax]:.4f}, "
          f"stress={s[imax]:.2f} GPa; "
          f"valley at strain={strain_eng[imin]:.4f}, "
          f"stress={s[imin]:.2f} GPa; drop={s[imax]-s[imin]:.2f} GPa")


# ========== 图 1：工程应力–应变 ==========
plt.figure(figsize=(10, 6))

# 原始数据（淡一点）
plt.plot(strain_eng, stress_eng, color="#bbbbbb", linewidth=1, alpha=0.7, label="Raw")

# 平滑数据（主曲线）
plt.plot(strain_eng, stress_eng_smooth, color="#1f77b4", linewidth=2.2, label="Smoothed")

# 拟合的弹性直线
fit_line = np.polyval(coeffs, strain_eng[mask_elastic])
plt.plot(strain_eng[mask_elastic], fit_line, "--", color="#ff7f0e", linewidth=2,
         label=f"Linear fit (E ≈ {E_fit:.1f} GPa)")

# 屈服点标记
plt.scatter(yield_strain, yield_stress, color="red", s=60, zorder=5, label="Yield point")
plt.annotate(f"Yield\n({yield_strain:.2f}, {yield_stress:.1f} GPa)",
             xy=(yield_strain, yield_stress),
             xytext=(yield_strain + 0.01, yield_stress + 0.5),
             arrowprops=dict(arrowstyle="->", color="red"),
             fontsize=10)

# 锯齿“谷值”位置标记
if len(serration_min_idx) > 0:
    serration_strains = strain_eng[serration_min_idx]
    serration_stresses = stress_eng_smooth[serration_min_idx]
    plt.scatter(serration_strains, serration_stresses, color="darkgreen",
                s=40, marker="v", zorder=5, label="Serration events")

plt.xlabel("Engineering strain", fontsize=14)
plt.ylabel("Engineering stress (GPa)", fontsize=14)
plt.title("Cu Engineering Stress–Strain Curve", fontsize=18, weight='bold')

plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.legend(fontsize=10)
plt.tight_layout()

plt.savefig("stress_strain_engineering_advanced.png", dpi=300)
plt.show()


# ========== 图 2：真应力–真应变 ==========
plt.figure(figsize=(10, 6))

# 原始真应力
plt.plot(strain_true, stress_true, color="#bbbbbb", linewidth=1, alpha=0.7, label="Raw true")

# 平滑后的真应力
plt.plot(strain_true, stress_true_smooth, color="#2ca02c", linewidth=2.2, label="Smoothed true")

# 把屈服点映射到真应力–真应变上
yield_strain_true = np.log(1.0 + yield_strain)
yield_stress_true = yield_stress * (1.0 + yield_strain)
plt.scatter(yield_strain_true, yield_stress_true, color="red", s=60, zorder=5, label="Yield point")

# 对应的 serration 位置在真应力–真应变图上（同样用谷值索引）
if len(serration_min_idx) > 0:
    serration_strain_true = np.log(1.0 + strain_eng[serration_min_idx])
    serration_stress_true = stress_true_smooth[serration_min_idx]
    plt.scatter(serration_strain_true, serration_stress_true, color="darkgreen",
                s=40, marker="v", zorder=5, label="Serration events")

plt.xlabel("True strain", fontsize=14)
plt.ylabel("True stress (GPa)", fontsize=14)
plt.title("Cu True Stress–Strain Curve", fontsize=18, weight='bold')

plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.legend(fontsize=10)
plt.tight_layout()

plt.savefig("stress_strain_true_advanced.png", dpi=300)
plt.show()

