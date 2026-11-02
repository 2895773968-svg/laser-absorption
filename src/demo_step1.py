# -*- coding: utf-8 -*-
"""
demo_step1.py —— 毕设第一步演示：用菲涅尔方程计算金属对激光的吸收率
=====================================================================
物理目标：
    给定金属的复折射率 N = n + i*k（n=折射率实部, k=消光系数），
    计算该金属表面对任意波长、任意入射角、任意偏振态激光的吸收率 A。
    金属在可见-近红外波段不透明（透射率 T≈0），故能量守恒简化为 A = 1 - R。

数据来源：
    refractiveindex.info 公开数据库（CC0 公有领域），
    Cu: P. B. Johnson & R. W. Christy, Phys. Rev. B 6, 4370 (1972)
    Al: A. D. Rakic, Appl. Opt. 34, 4755 (1995)
    Fe: P. B. Johnson & R. W. Christy, Phys. Rev. B 9, 5056 (1974)
"""

import numpy as np                  # 向量化数值计算核心库
import yaml                         # 解析 refractiveindex.info 的 YAML 数据文件
import matplotlib.pyplot as plt     # 绘图
import os

# ---------- 让 matplotlib 正常显示中文（Windows 自带微软雅黑） ----------
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False   # 修复负号显示为方块的字体问题


# =====================================================================
# 第一部分：数据读取 —— 把 YAML 文本变成 numpy 数组
# =====================================================================
def load_nk_yaml(path):
    """读取 refractiveindex.info 格式的 nk 数据文件。

    YAML 结构（已用真实文件验证）：
        DATA:
          - type: tabulated nk
            data: |
                <波长(微米)> <n> <k>     <- 每行三个数，空格分隔
    返回：
        wl : np.ndarray  波长数组，单位换算为纳米 (nm)
        n  : np.ndarray  折射率实部
        k  : np.ndarray  消光系数（虚部量值）
    """
    with open(path, "r", encoding="utf-8") as f:
        doc = yaml.safe_load(f)              # 解析整个 YAML 为 Python 字典

    # 数据库文件中 DATA 是一个列表，通常第一个元素就是 tabulated nk 表
    raw = doc["DATA"][0]["data"]             # 取出原始多行字符串
    rows = np.array([
        [float(x) for x in line.split()]     # 每行按空格切成 3 个浮点数
        for line in raw.strip().splitlines() # strip() 去掉首尾空行
    ])

    # 第 0 列是波长，单位微米 -> 乘 1000 转纳米，与"激光波长"习惯一致
    wl = rows[:, 0] * 1000.0
    n, k = rows[:, 1], rows[:, 2]            # 第 1、2 列分别是 n 和 k

    # 保证波长严格递增（个别文件可能有乱序），否则 np.interp 会出错
    order = np.argsort(wl)
    return wl[order], n[order], k[order]


def get_nk(wl_data, n_data, k_data, wl_query):
    """在测量点上做线性插值，得到任意波长处的 (n, k)。

    物理注解：数据库只在离散波长点给出实验值；
    金属光学常数在无带间跃迁的平滑区随波长缓变，
    线性插值在毕设精度下完全够用（更严格可做样条或 Drude 拟合）。
    """
    # np.clip 把查询点限制在数据范围内，避免外推越界
    wl_q = np.clip(np.asarray(wl_query, dtype=float),
                   wl_data.min(), wl_data.max())
    return (np.interp(wl_q, wl_data, n_data),   # n(λ) 线性插值
            np.interp(wl_q, wl_data, k_data))   # k(λ) 线性插值


# =====================================================================
# 第二部分：物理计算 —— 菲涅尔方程（复折射率版）
# =====================================================================
def normal_absorptance(n, k):
    """垂直入射（θ=0）时的吸收率解析式，可直接验证程序正确性。

    推导：垂直入射时振幅反射系数 r = (1 - N)/(1 + N)，
          R = |r|^2 = [(n-1)^2 + k^2] / [(n+1)^2 + k^2]，
          A = 1 - R = 4n / [(n+1)^2 + k^2]。
    注意 k 出现在分母的 k^2 项——这就是为什么 Cu 的 k 在红外高达 ~6.5，
    吸收率反而被"压"到只有百分之几：k 越大，电子再辐射（反射）越强。
    """
    return 4.0 * n / ((n + 1.0) ** 2 + k ** 2)


def fresnel_absorptance(n, k, theta_deg, polarization="unpolarized"):
    """斜入射菲涅尔方程（吸收介质版），支持向量化。

    物理模型（电磁场边界条件推导，见玻恩《光学原理》§1.6 / §13）：
      入射介质为空气（近似真空, N1 = 1），金属复折射率 N2 = n + i*k。
      复介电函数 eps = N2^2 = (n^2 - k^2) + i*2nk。
      斜入射时折射角为复数，为避免显式处理它，引入中间量：
          u = sqrt(eps - sin^2(theta))
      它等价于 N2*cos(theta_t)（复折射角余弦），是标准有效形式。

    参数:
        theta_deg    : 入射角（度），可以是标量或 numpy 数组
        polarization : 's' / 'p' / 'unpolarized'（非偏振取 (A_s+A_p)/2）
    返回:
        A : 吸收率数组（金属不透明，T≈0，故 A = 1 - R）
    """
    theta = np.deg2rad(np.asarray(theta_deg, dtype=float))

    eps = (n + 1j * k) ** 2                          # 复介电函数 ε = N²
    sin2 = np.sin(theta) ** 2                        # sin²θ，实数
    u = np.sqrt(eps - sin2)                          # 关键中间量（复数开方）

    cos_t = np.cos(theta)                            # 入射侧方向余弦 cosθ（实数）

    # --- 菲涅尔振幅反射系数（有效形式，与教科书 r_s/r_p 完全等价） ---
    r_s = (cos_t - u) / (cos_t + u)                  # s 偏振（E⊥入射面）
    r_p = (eps * cos_t - u) / (eps * cos_t + u)      # p 偏振（E∥入射面）

    # R = |r|²：r 是复数，模方即功率反射率
    A_s = 1.0 - np.abs(r_s) ** 2
    A_p = 1.0 - np.abs(r_p) ** 2

    if polarization == "s":
        return A_s
    if polarization == "p":
        return A_p
    return 0.5 * (A_s + A_p)                         # 非偏振光 = 两个偏振平均


# =====================================================================
# 第三部分：演示入口 —— 计算三种金属 + 出两张图
# =====================================================================
if __name__ == "__main__":
    base = os.path.join(os.path.dirname(__file__), "..", "data")
    out = os.path.join(os.path.dirname(__file__), "..", "output")
    os.makedirs(out, exist_ok=True)

    # ---- 1. 载入三种金属的光学常数 ----
    metals = {
        "Cu": load_nk_yaml(os.path.join(base, "Cu_Johnson.yml")),
        "Al": load_nk_yaml(os.path.join(base, "Al_Rakic.yml")),
        "Fe": load_nk_yaml(os.path.join(base, "Fe_Johnson.yml")),
    }

    # ---- 2. 工业常用激光波长处的吸收率表（垂直入射） ----
    laser_lines = [450, 515, 532, 808, 980, 1064]   # 蓝/绿/红/近红外常见激光
    print("=" * 64)
    print("垂直入射吸收率 A(%) —— 正入射解析式 A = 4n/[(n+1)²+k²]")
    print(f"{'波长(nm)':>10} | {'Cu':>8} | {'Al':>8} | {'Fe':>8}")
    print("-" * 64)
    for wl in laser_lines:
        row = []
        for m, (wd, nd, kd) in metals.items():
            n_q, k_q = get_nk(wd, nd, kd, wl)
            row.append(normal_absorptance(n_q, k_q) * 100)
        print(f"{wl:>10} | {row[0]:>8.2f} | {row[1]:>8.2f} | {row[2]:>8.2f}")
    print("=" * 64)
    print("预期对照：Cu 在 1064nm 约 2~5%，在 532nm 跳到 40~60%，"
          "在 450nm 达 60%+（d 带带间跃迁）—— 若对上说明物理正确")

    # ---- 3. 图 1：吸收率 vs 波长（400-1100 nm 扫描） ----
    lam = np.linspace(400, 1100, 701)                # 1 nm 步长扫描
    fig, ax = plt.subplots(figsize=(8, 5), dpi=110)
    for (name, (wd, nd, kd)), color in zip(metals.items(),
                                           ["#BA7517", "#888780", "#185FA5"]):
        n_q, k_q = get_nk(wd, nd, kd, lam)
        ax.plot(lam, normal_absorptance(n_q, k_q) * 100,
                label=f"{name}", color=color, lw=2)
    # 标出蓝/绿/红外三条工业激光线，直观呼应上一轮讨论
    for lw, c, txt in [(450, "#378ADD", "蓝光450"), (532, "#639922", "绿光532"),
                       (1064, "#A32D2D", "光纤1064")]:
        ax.axvline(lw, ls="--", lw=0.8, color=c, alpha=0.6)
        ax.text(lw + 6, 68, txt, fontsize=9, color=c)
    ax.set_xlabel("波长 λ (nm)")
    ax.set_ylabel("吸收率 A (%)")
    ax.set_title("金属吸收率随波长变化 —— 垂直入射，菲涅尔方程计算\n"
                 "数据: refractiveindex.info (Johnson-Christy 1972/74, Rakic 1995)")
    ax.legend(title="金属")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "step1_A_vs_wavelength.png"))
    print(f"已保存: {os.path.join(out, 'step1_A_vs_wavelength.png')}")

    # ---- 4. 图 2：斜入射效应 —— Cu@1064nm 的 s/p 偏振吸收率 vs 入射角 ----
    # 物理预期：p 偏振在 60-80° 附近出现类似"布儒斯特角"的吸收峰，
    # 这是斜入射切割/焊接中" graze 角耦合增强"现象的来源
    n_cu, k_cu = get_nk(*metals["Cu"], 1064.0)       # Cu 在 1064nm 处的 n,k
    angles = np.linspace(0, 89, 178)                  # 0~89 度扫描
    A_s = fresnel_absorptance(n_cu, k_cu, angles, "s") * 100
    A_p = fresnel_absorptance(n_cu, k_cu, angles, "p") * 100
    A_un = fresnel_absorptance(n_cu, k_cu, angles, "unpolarized") * 100

    fig2, ax2 = plt.subplots(figsize=(8, 5), dpi=110)
    ax2.plot(angles, A_s, label="s 偏振", color="#185FA5", lw=2)
    ax2.plot(angles, A_p, label="p 偏振", color="#D85A30", lw=2)
    ax2.plot(angles, A_un, label="非偏振(平均)", color="#444441", lw=2, ls="--")
    ax2.set_xlabel("入射角 θ (°)")
    ax2.set_ylabel("吸收率 A (%)")
    ax2.set_title(f"Cu 在 1064 nm 的斜入射吸收（n={n_cu:.3f}, k={k_cu:.3f}）\n"
                  "p 偏振吸收峰 = 金属版'布儒斯特效应'")
    ax2.legend()
    ax2.grid(alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out, "step1_Cu_angle_1064.png"))
    print(f"已保存: {os.path.join(out, 'step1_Cu_angle_1064.png')}")
