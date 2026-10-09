# H桥单极性SPWM 变频器输出谐波分析 (MATLAB)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2025b%20(R2019b%2B)-orange.svg)](#7-环境)

把原来的 Python 发波仿真脚本 `h_bridge_svpwm.py` 移植到 MATLAB，并加上**精确的谐波成分计算**。
调制规律与原 Python 脚本**逐点一致**（三角载波、单极性SPWM、H桥输出 `Vout = Vdc*(Q1-Q3)`），
但谐波不是靠"1 µs 采样 + FFT"估算，而是**解析求解每个开关跳变时刻后逐段精确积分**，
结果精确到机器精度（1e-13 量级），可以作为标准答案。

---

## 1. 文件清单

| 文件 | 作用 |
|---|---|
| `run_harmonic_analysis.m` | **主脚本**：改参数 → 按 F5 运行，出表格 + 图 + CSV |
| `hbridge_spwm_harmonics.m` | **核心函数**：一次算完某工况的全部谐波（精确法 + 采样FFT法） |
| `verify_harmonics.m` | **自检脚本**：13 项验证，含 4 种独立方法的交叉验证 |
| `spwm_tri_carrier.m` / `spwm_out_level.m` / `spwm_level_residual.m` | 载波与电平函数（与 Python 逐点一致） |
| `harmonic_report.m` | 打印谐波报表（低次谐波表、最大谐波、载波谐波群） |
| `sweep_modulation_index.m` / `plot_thd_vs_M.m` | 调制比 M = 0.1~1.0 扫描与作图 |
| `plot_waveforms.m` | 图1：载波/调制波、4路驱动、输出电压（含开关细节） |
| `plot_spectrum.m` | 图2：幅值谱（线性/对数）+ 与原 1 µs 采样 FFT 的误差 |
| `plot_spectrum_detail.m` | 图3：0~2 kHz、2fc/4fc 谐波群、最大 20 个谐波 |
| `plot_spectra_vs_M.m` | 图5：不同 M 下的频谱对比 |
| `export_spectrum_csv.m` / `export_fig_png.m` | 导出 CSV / PNG |
| `output/` | 运行后生成的数据与图形 |

---

## 2. 快速使用

```matlab
% 1) 先跑一次自检（建议第一次使用时执行）
verify_harmonics          % 13 项全部 PASS 才算环境正常

% 2) 全流程分析（改 run_harmonic_analysis.m 顶部的参数后按 F5）
run_harmonic_analysis     % 生成 output/ 下的 5 张图 + 2 个 CSV

% 3) 单独算某个工况（返回结构体 R，字段见函数头注释）
R = hbridge_spwm_harmonics(2400, 50, 100, 0.8);
R.Amp(95)                 % 4750 Hz 谐波幅值
R.THD_all, R.WTHD         % 畸变指标
plot_spectrum(R)          % 画频谱
```

### 主要参数（`run_harmonic_analysis.m` 顶部）

| 参数 | 默认 | 说明 |
|---|---|---|
| `fc` | 2400 | 开关/载波频率 Hz（原 Python 是 400，可调；**建议取 f0 的整数倍**） |
| `f0` | 50 | 调制波（基波）频率 Hz |
| `Vdc` | 100 | 直流母线电压 V |
| `M` | 0.8 | 调制比 0.1 ~ 1.0（线性调制区） |
| `Ncyc` | 10 | 参与分析的基波周期数 |
| `DT` | 1e-6 | 采样步长（仅用于复现原 Python 的 FFT 对比） |
| `Kmax` | 4000 | 最高谐波次数，4000×50 = 200 kHz |
| `Rload` / `Lload` | 0.5 / 2e-3 | 负载 R-L，用于算**电流**谐波（设 `Rload=0` 则只算电压） |

---

## 3. 两个计算方法

1. **精确解析法（推荐，本工具的主要结果）**
   - 输出电平只可能是 `+Vdc / 0 / -Vdc`（`|v|>|c|` 时输出 `sign(v)*Vdc`，否则为 0，即续流）。
   - 令 `F(t) = |v(t)| - |c(t)| = 0`，用向量化二分法求出**每个跳变时刻**（收敛到机器精度）。
   - 把波形写成"分段常数"，每段解析积分得到傅里叶系数：
     `c_h = (1/T)·Σ_j lev_j·(b_j-a_j)·e^{-jπ h f0 (a_j+b_j)}·sinc(h f0 (b_j-a_j))`
   - 与时间网格无关，**没有量化误差**。

2. **采样 FFT 法（复现原 Python）**
   - 完全照原脚本 `DT = 1 µs` 等步长采样后 FFT，用来评估**原仿真的精度损失**。
   - 采样点不可能正好落在开关时刻上 → 每个跳变时刻有最多 ±1 µs 的抖动。

---

## 4. 结果与物理结论（fc = 2400 Hz, f0 = 50 Hz, Vdc = 100 V, M = 0.8）

### 4.1 基波
**A₁ = M·Vdc = 80.000000 V**（计算误差 1.4e-14 V，与理论完全一致）。
谐波分析里 A₁ 随 M 严格线性，见 `fig4` 子图 (a)。

### 4.2 谐波分布规律（单极性SPWM的固有特性）

* **只有奇次谐波**（波形满足半波对称，fc/f0 为整数时成立）。
* **低于 2fc 的谐波全部严格为 0**：第 3、5、7、9… 次谐波幅值 < 2e-12 V，即"理论上不存在"。
  → 这是单极性SPWM最大的优点：**不需要很大的输出滤波器就能得到低失真的正弦电流**。
* 谐波只出现在**偶数倍载波**附近：`h = 96m ± (2k-1)`（96 = 2fc/f0），边带间隔 = 2f₀ = 100 Hz，
  中心点 h = 96（正好 4800 Hz）为 0。
* 理论幅值有闭式解（本工具已与它对比验证到 1e-11）：

  ```
  A(m,n) = Vdc · (4/π) · (1/m) · |J_n(m·π·M/2)|      m = 2,4,6,…   n = ±1,±3,…
  ```

**M = 0.8 时最主要的谐波（峰值 V / 占基波 %）：**

| 频率 | h | 幅值 | 占基波 |
|---|---|---|---|
| 4750 Hz (2fc−f0) | 95 | 31.435 | 39.29 % |
| 4850 Hz (2fc+f0) | 97 | 31.435 | 39.29 % |
| 4650 / 4950 Hz | 93 / 99 | 13.947 | 17.43 % |
| 9450 / 9750 Hz (4fc±) | 189 / 195 | 11.465 | 14.33 % |
| 9550 / 9650 Hz | 191 / 193 | 10.518 | 13.15 % |
| 4800 Hz (2fc 中心) | 96 | 0 | 0 % |

### 4.3 畸变指标

⚠️ **电压 THD 强烈依赖求和上限**，报告数字时必须说明带宽，否则没有意义。本工具同时给出：

| 求和范围 | M=0.8 |
|---|---|
| h = 2~50（2.5 kHz） | **0.000 %** |
| h ≤ 2fc（4.8 kHz） | 43.017 % |
| h ≤ 4fc（9.6 kHz） | 64.768 % |
| h ≤ 200 kHz | 76.433 % |
| WTHD（1/h 加权，更能反映电流畸变） | **0.660 %** |
| 负载电流 THD（R=0.5 Ω, L=2 mH） | 0.843 % |

**调制比越小，电压 THD 越大**（因为基波 A₁=M·Vdc 变小，而开关谐波幅值基本不变）：

| M | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| THD ≤2fc / % | 98.8 | 95.2 | 89.4 | 81.7 | 72.7 | 62.8 | 52.7 | 43.0 | 34.6 | 28.1 |
| THD ≤100 kHz / % | 327.9 | 226.3 | 177.1 | 145.7 | 122.8 | 104.6 | 89.4 | 75.9 | 63.5 | 51.4 |
| WTHD / % | 1.73 | 1.57 | 1.41 | 1.26 | 1.10 | 0.95 | 0.80 | 0.66 | 0.53 | 0.43 |
| 电流 THD ≤100 kHz / % | 2.21 | 2.01 | 1.80 | 1.60 | 1.41 | 1.21 | 1.02 | 0.84 | 0.68 | 0.55 |

> 结论：**低 M 时不要看电压 THD**（那是"基波太小"造成的假象），要看 WTHD 或电感负载的电流 THD。

### 4.4 原 Python 1 µs 仿真的精度警告

在 fc = 2400 Hz 下，用 `DT = 1 µs` 等步长采样做 FFT（即原 Python 的做法），实测误差：

| 谐波 | 精确解 | 1 µs 采样 FFT | 误差 |
|---|---|---|---|
| 50 Hz（基波） | 80.0000 V | 79.9332 V | **−0.084 %** |
| 150 Hz（h=3） | ~0（1e-13） | **0.0327 V** | 虚假谐波 |
| 250 Hz（h=5） | ~0 | 0.0295 V | 虚假谐波 |
| 350 Hz（h=7） | ~0 | 0.0453 V | 虚假谐波 |
| 4750 Hz（h=95） | 31.4353 V | 31.5528 V | +0.374 % |
| 4850 Hz（h=97） | 31.4353 V | 31.3643 V | −0.226 % |
| 显著谐波（>1 % Vdc）最大 | — | — | **10.2 %**（h=471, 23.55 kHz） |
| 最大绝对误差 | — | — | 0.189 V（0.19 % Vdc） |

**两个要特别注意的后果：**

1. 用 1 µs 仿真去验证"低次谐波为零"会看到一个 **0.01~0.05 V 的虚假本底**
   （约 0.01~0.06 % 基波），而真实值应该是 1e-13 V 量级。别把它当成真实的低次谐波。
2. 2fc 谐波群（h = 91…101）显著谐波的幅值误差最高约 **3.9 %**，用来看"开关频率谐波有多大"
   勉强可以，用来对比不同调制策略/滤波器效果就不可靠了。

原因：1 µs 采样使每个开关沿抖动最多 ±1 µs，而载波周期只有 417 µs（抖动 0.24 %），
对 4750 Hz 谐波相当于 0.03 rad 的相位误差。
**建议：原 Python 脚本若要算频谱，把 `DT` 减到 1e-8 以下，或直接用本工具的精确法。**

---

## 5. 验证（`verify_harmonics.m`，13 项全部 PASS）

| 项 | 检查内容 | 实测 |
|---|---|---|
| T1 | 基波 = M·Vdc | 1.99e-13 V |
| T2 | 直流分量为 0 | 1.23e-13 V |
| T3 | 偶次谐波为 0 | 1.56e-12 V |
| T4 | 3~95 次奇次谐波为 0 | 1.89e-12 V |
| T5 | 相位 = ±90°（奇函数，余弦分量=0） | 7.9e-12 V |
| T6 | 帕塞瓦尔定理（波形/频谱有效值） | 截断差 0.24 % |
| T7 | 1 µs 采样 FFT 的有效值 | 误差 0.058 % |
| T8 | vs **10 ns 细采样 FFT**（独立数值方法） | 0.21 % |
| T9 | vs **16 点 Gauss-Legendre 数值积分**（独立数值方法） | 4.2e-14 |
| T10 | A₁ 随 M 严格线性（斜率 = Vdc） | 5.6e-13 V |
| T11 | 奇数载波比 fc/f0 = 49 仍正确 | 通过 |
| T12 | M = 0.02 / 0.1 / 1.0 边界 | 通过 |
| T13 | vs **双重傅里叶级数（贝塞尔函数）解析解**（独立理论，94 个谐波） | 1.1e-11 |

T13 是最强的一项：它用完全不同的数学路径（Bessel 函数解析式）复核了所有载波谐波群的幅值。

---

## 6. 常见问题

**Q: 改 fc 后结果不对 / 频谱出现"毛刺"？**
要求 `fc/f0` 为整数（如 2400/50 = 48）。若不为整数，程序会给出警告——此时波形在 1/f0 内
不是严格周期的，任何 FFT 都会泄漏。可令 `Ncyc` 使 `Ncyc*fc/f0` 为整数，或直接取 fc 为 f0 的整数倍。

**Q: 想算三相变频器的线电压谐波？**
把 `spwm_out_level.m` 换成三相调制（SVPWM/SPWM + 第三谐波注入）并取线电压 `Vab = Va - Vb`，
其余流程（精确开关时刻 + 解析积分）完全可复用。

**Q: 想算电机电流的谐波？**
设 `Rload`、`Lload` 为电机等效参数，程序会按 `I_h = V_h / |R + j·2πh·f0·L|` 给出电流各次谐波
及电流 THD（更贴近电机的实际发热与转矩脉动）。也可以直接对定位的电流波形做 FFT。

**Q: 为什么和示波器/实测的频谱不完全一致？**
本工具是**理想开关**模型（与原 Python 一致）：没有死区、没有开关上升/下降时间、没有管压降、
没有母线波动、没有共模与寄生参数。实机上：
* **死区**会引入 3、5、7 次低次谐波（本工具中它们严格为 0）——这是实测与仿真最大的差异来源；
* 开关沿的 dv/dt 有限会削弱高次谐波（> 50 次以后）；
* 负载/电缆的寄生电容电感会形成谐振尖峰。
如果两者差异很大，优先检查死区设置与开关沿时间。

**Q: 输出文件分别是什么？**
* `spectrum_fc2400_M0.80.csv`：h、频率、精确幅值、占基波%、相位、1 µs FFT 幅值及其误差、电流幅值；
* `thd_vs_M.csv`：调制比扫描表；
* `fig1`~`fig5`：波形、频谱、频谱细节、THD-M 曲线、不同 M 的频谱对比。

---

## 7. 环境

MATLAB R2025b（R2019b 以上应该都可以，用到 `arguments` 名称-值语法、`exportgraphics`、`tiledlayout`）。
不需要任何工具箱（`besselj` 属于基础 MATLAB）。

---

## 8. 许可证 (License)

本项目以 **GNU General Public License v3.0 或更新版本（GPL-3.0-or-later）** 发布，完整条款见 [LICENSE](LICENSE)。

```
hbridge-spwm-harmonics —— H桥单极性SPWM输出电压谐波精确分析
Copyright (C) 2026 leuk0cyte

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
```

每个源文件头部均带 SPDX 标识：

```matlab
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
```

### 贡献与免责

* 欢迎提交 Issue / Pull Request；按照 GPLv3 第 5 条，修改后的版本需同样以 GPLv3 发布并注明修改。
* 本工具是**理想开关**的仿真与解析计算程序，用于谐波分析研究与教学，**不构成任何工程担保**。
  用于实际变频器设计（死区、开关损耗、母线波动、EMC）时请自行验证，见 [第 6 节](#6-常见问题)。

