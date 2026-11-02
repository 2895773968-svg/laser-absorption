# -*- coding: utf-8 -*-
"""
material_optics.py —— MaterialOptics 类
========================================
单一金属的光学常数容器 + 菲涅尔吸收率计算器。

职责边界（单一职责原则）：
    只负责 "读取数据(n,k) + 物理计算(吸收率)"，不涉及任何绘图。
    绘图统一交给 Plotter 类，保证物理内核可独立单元测试。

物理依据：
    复折射率 N = n + i*k，复介电函数 eps = N^2；
    斜入射菲涅尔方程（有效形式）：
        u   = sqrt(eps - sin^2(theta))          # = N2*cos(theta_t)，复数
        r_s = (cos(theta) - u) / (cos(theta) + u)
        r_p = (eps*cos(theta) - u) / (eps*cos(theta) + u)
        A   = 1 - |r|^2   （金属不透明，透射率 T ≈ 0）
"""

import numpy as np
import yaml


# =====================================================================
# 纯函数版菲涅尔计算 —— 放在类外，方便 pytest 直接对拍教科书公式
# =====================================================================
def fresnel_absorptance(n, k, theta_deg, polarization="unpolarized"):
    """斜入射菲涅尔吸收率，支持 numpy 广播。

    参数
    ----
    n, k         : 复折射率实部/消光系数。可以是标量、一维数组（随波长）。
    theta_deg    : 入射角（度）。可以是标量或一维数组（随角度）。
                   当 n 和 theta 都是一维数组时，返回形状 (len(n), len(theta)) 的二维网格；
                   其余混合情形由 numpy 广播规则自然处理。
    polarization : 's' | 'p' | 'unpolarized'

    返回
    ----
    A : 吸收率数组（无量纲，0~1）
    """
    theta = np.deg2rad(np.atleast_1d(np.asarray(theta_deg, dtype=float)))
    eps = (np.atleast_1d(n) + 1j * np.atleast_1d(k)) ** 2   # 复介电函数
    sin2 = np.sin(theta) ** 2                                # sin²θ，实数

    # 二维网格情形：波长轴放在第 0 维、角度轴放在第 1 维
    grid = (eps.size > 1) and (sin2.size > 1)
    if grid:
        e, s = eps[:, None], sin2[None, :]          # (Nλ,1) 与 (1,Nθ)
        cos_t = np.cos(theta)[None, :]              # cosθ 广播到 (1,Nθ)
    else:
        e, s, cos_t = eps, sin2, np.cos(theta)      # 一维混合广播

    u = np.sqrt(e - s)                              # 关键中间量（复数开方）

    r_s = (cos_t - u) / (cos_t + u)                 # s 偏振振幅反射系数
    r_p = (e * cos_t - u) / (e * cos_t + u)         # p 偏振振幅反射系数

    A_s = 1.0 - np.abs(r_s) ** 2                    # 功率反射率 R=|r|²，A=1-R
    A_p = 1.0 - np.abs(r_p) ** 2

    if polarization == "s":
        return A_s
    if polarization == "p":
        return A_p
    if polarization == "unpolarized":
        return 0.5 * (A_s + A_p)                    # 非偏振 = 两种偏振平均
    raise ValueError(f"未知偏振类型: {polarization!r}，应为 's'/'p'/'unpolarized'")


# =====================================================================
# MaterialOptics 主类
# =====================================================================
class MaterialOptics:
    """一种金属的 n(λ), k(λ) 数据容器 + 吸收率计算接口。"""

    def __init__(self, name, wavelength_nm, n, k, reference=""):
        """
        参数
        ----
        name         : 金属符号，如 'Cu'
        wavelength_nm: 波长数组（nm，严格递增）
        n, k         : 与 wavelength_nm 等长的光学常数数组
        reference    : 数据来源描述（写入图表标注 / README）
        """
        self.name = name
        self.reference = reference
        order = np.argsort(wavelength_nm)            # 保证波长递增（插值前提）
        self.wl = np.asarray(wavelength_nm, float)[order]
        self.n = np.asarray(n, float)[order]
        self.k = np.asarray(k, float)[order]

    # ---------- 构造方法 ----------
    @classmethod
    def from_yaml(cls, name, path, reference=""):
        """从 refractiveindex.info 的 YAML 文件构建实例。

        YAML 核心结构：
            DATA:
              - type: tabulated nk
                data: |
                    <波长(微米)> <n> <k>    # 每行三列，空格分隔
        """
        with open(path, "r", encoding="utf-8") as f:
            doc = yaml.safe_load(f)
        rows = np.array([[float(x) for x in line.split()]
                         for line in doc["DATA"][0]["data"].strip().splitlines()])
        if not reference:
            # 自动把 YAML 里的 REFERENCES 字段压缩成一行，便于溯源
            reference = " ".join(str(doc.get("REFERENCES", "")).split())[:120]
        return cls(name, rows[:, 0] * 1000.0,       # 微米 -> 纳米
                   rows[:, 1], rows[:, 2], reference)

    # ---------- 数据接口 ----------
    def nk(self, wl):
        """任意波长处的 (n, k)，线性插值（带越界截断保护）。"""
        wl_q = np.clip(np.asarray(wl, float), self.wl.min(), self.wl.max())
        return np.interp(wl_q, self.wl, self.n), np.interp(wl_q, self.wl, self.k)

    # ---------- 物理计算接口 ----------
    def absorptance_normal(self, wl):
        """垂直入射吸收率解析式 A = 4n/[(n+1)²+k²]（用于快速计算与代码自检）。"""
        n, k = self.nk(wl)
        return 4.0 * n / ((n + 1.0) ** 2 + k ** 2)

    def absorptance(self, wl, theta_deg=0.0, polarization="unpolarized"):
        """通用吸收率：波长与角度均可数组化，返回广播后的数组。"""
        n, k = self.nk(wl)
        return fresnel_absorptance(n, k, theta_deg, polarization)

    def __repr__(self):
        return (f"MaterialOptics({self.name}, {self.wl.min():.0f}-"
                f"{self.wl.max():.0f} nm, {self.wl.size} points)")
