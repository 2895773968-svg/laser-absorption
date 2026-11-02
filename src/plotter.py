# -*- coding: utf-8 -*-
"""
plotter.py —— Plotter 类
=========================
学术发表级图表生成器。

风格约定（publication style）：
    - 衬线字体 Times New Roman + STIX 数学字体（物理期刊主流风格）
    - 300 DPI，PNG（README 用）+ PDF（矢量，论文投稿用）双格式输出
    - 轴标签带物理量与单位：Wavelength (nm) / Incident angle (°) / A (%)
    - 色板遵循"暖色=强吸收（热）"直觉：Cu 琥珀、Al 灰、Fe 蓝
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from material_optics import MaterialOptics
from angle_sweep import AngleSweepResult

# 每种金属的固定配色（保证全项目图表颜色一致，可复现）
METAL_COLORS = {"Cu": "#BA7517", "Al": "#888780", "Fe": "#185FA5",
                "Au": "#B8860B", "Ag": "#7F7F7F"}


class Plotter:
    """消费 MaterialOptics / AngleSweepResult，产出发表级图表。"""

    # ---- 发表级 rcParams（一次设置，全局生效） ----
    PUBLICATION_RC = {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "STIXGeneral"],
        "mathtext.fontset": "stix",          # 数学字体与正文匹配（物理期刊风格）
        "axes.unicode_minus": False,
        "axes.linewidth": 0.8,               # 细轴线，期刊惯例
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "legend.frameon": False,             # 图例无边框，更干净
        "figure.dpi": 110,
        "savefig.dpi": 300,                  # 保存即 300 DPI，满足期刊/印刷要求
        "savefig.bbox": "tight",
    }

    def __init__(self, style="publication", out_dir="output"):
        self.rc = dict(self.PUBLICATION_RC) if style == "publication" else {}
        self.out_dir = out_dir
        os.makedirs(out_dir, exist_ok=True)

    def _finalize(self, fig, filename):
        """统一保存：PNG(300DPI, README 用) + PDF(矢量, 论文用)。"""
        for ext in (".png", ".pdf"):
            fig.savefig(os.path.join(self.out_dir, filename + ext))
        plt.show() if plt.isinteractive() else plt.close(fig)
        return os.path.join(self.out_dir, filename + ".png")

    # ------------------------------------------------------------------
    def plot_spectrum(self, materials, wl_range=(400, 1100), num=701,
                      theta_deg=0.0, polarization="unpolarized",
                      laser_lines=(450, 532, 1064), filename="spectrum"):
        """图 A：吸收率 vs 波长（可叠加多种金属对比）。

        参数
        ----
        materials : list[MaterialOptics]，至少一种金属
        """
        with plt.rc_context(self.rc):
            lam = np.linspace(*wl_range, num)
            fig, ax = plt.subplots(figsize=(5.2, 3.6))
            for m in materials:
                A = m.absorptance(lam, theta_deg, polarization) * 100
                ax.plot(lam, A, lw=1.8, color=METAL_COLORS.get(m.name, None),
                        label=m.name)
            # 标注工业激光波长（虚线竖线 + 顶部文字）
            for lw_, txt in zip(laser_lines,
                                [f"{w} nm" for w in laser_lines]):
                ax.axvline(lw_, ls="--", lw=0.7, color="grey", alpha=0.6)
                ax.text(lw_ + 6, 0.97, txt, transform=ax.get_xaxis_transform(),
                        fontsize=8, color="grey")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel(f"Absorptance A (%)")       # 物理量 + 单位
            pol_txt = {"s": "s-pol", "p": "p-pol",
                       "unpolarized": "unpolarized"}[polarization]
            ax.set_title(f"Fresnel absorptance, {pol_txt}, "
                         f"θ = {theta_deg:.0f}°")
            ax.legend(title="Metal")
            ax.grid(alpha=0.3, lw=0.5)
            fig.tight_layout()
            return self._finalize(fig, filename)

    # ------------------------------------------------------------------
    def plot_angle_sweep(self, result: AngleSweepResult, filename=None):
        """图 B：固定波长下 A(θ) —— s/p/非偏振三条曲线 + 伪布儒斯特峰标注。"""
        if filename is None:
            filename = f"angle_{result.material_name}_{result.wavelength_nm:.0f}"
        with plt.rc_context(self.rc):
            fig, ax = plt.subplots(figsize=(5.2, 3.6))
            ax.plot(result.angles_deg, result.A_s * 100, lw=1.8,
                    color="#185FA5", label="s-polarized")
            ax.plot(result.angles_deg, result.A_p * 100, lw=1.8,
                    color="#D85A30", label="p-polarized")
            ax.plot(result.angles_deg, result.A_un * 100, lw=1.4,
                    color="#444441", ls="--", label="unpolarized")
            # 标注 p 偏振峰值（伪布儒斯特角）
            A_pk, th_pk = result.p_peak
            ax.plot(th_pk, A_pk * 100, "v", color="#D85A30", ms=7)
            ax.annotate(f"pseudo-Brewster\n({th_pk:.1f}°, {A_pk*100:.1f}%)",
                        xy=(th_pk, A_pk * 100),
                        xytext=(th_pk - 28, A_pk * 100 + 3),
                        fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7))
            ax.set_xlabel("Incident angle (°)")
            ax.set_ylabel("Absorptance A (%)")
            ax.set_title(f"{result.material_name}, "
                         f"λ = {result.wavelength_nm:.0f} nm")
            ax.legend()
            ax.grid(alpha=0.3, lw=0.5)
            fig.tight_layout()
            return self._finalize(fig, filename)

    # ------------------------------------------------------------------
    def plot_heatmap(self, material: MaterialOptics,
                     wl_range=(400, 1100), angle_range=(0, 89),
                     n_wl=701, n_ang=179, polarization="p",
                     laser_lines=(450, 532, 1064),
                     cmap="inferno", filename=None):
        """图 C（README 主视觉）：波长 × 入射角二维吸收率热力图。

        实现：contourf 填充等值线热力图（比 imshow 更平滑、期刊更常用），
        自动标注全局最大值点（p 偏振伪布儒斯特吸收的波长-角度位置）。
        """
        if filename is None:
            filename = f"heatmap_{material.name}_{polarization}pol"
        # ---- 1. 构建二维物理网格并计算（MaterialOptics 自动广播） ----
        lam = np.linspace(*wl_range, n_wl)              # 波长轴 (Nλ,)
        ang = np.linspace(*angle_range, n_ang)          # 角度轴 (Nθ,)
        A = material.absorptance(lam, ang, polarization)  # (Nλ, Nθ)
        A_pct = A * 100

        # ---- 2. 找全局最大值点 ----
        i, j = np.unravel_index(np.argmax(A_pct), A_pct.shape)
        lam_pk, ang_pk, A_pk = lam[i], ang[j], A_pct[i, j]

        # ---- 3. 绘图（注意 Z 的形状须为 (Nθ, Nλ)，故转置） ----
        with plt.rc_context(self.rc):
            fig, ax = plt.subplots(figsize=(6.4, 4.2))
            levels = np.linspace(0, np.ceil(A_pct.max() / 10) * 10, 21)
            cf = ax.contourf(lam, ang, A_pct.T, levels=levels, cmap=cmap,
                             extend="max")
            cbar = fig.colorbar(cf, ax=ax, pad=0.02)
            cbar.set_label(f"Absorptance A (%)")

            # 激光线参考（白色虚线，工业波长定位）
            for w in laser_lines:
                ax.axvline(w, color="white", ls=":", lw=0.8, alpha=0.7)
                ax.text(w + 4, angle_range[1] * 0.04, f"{w}",
                        color="white", fontsize=8)

            # 最大值标注：星号 + 引出线
            ax.plot(lam_pk, ang_pk, "*", color="white", ms=13,
                    mec="black", mew=0.5)
            ax.annotate(f"max {A_pk:.1f}%\n({lam_pk:.0f} nm, {ang_pk:.1f}°)",
                        xy=(lam_pk, ang_pk),
                        xytext=(lam_pk + 60, ang_pk - 14),
                        color="white", fontsize=9,
                        arrowprops=dict(arrowstyle="->", color="white", lw=0.8))

            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Incident angle (°)")
            pol = {"s": "s-polarized", "p": "p-polarized",
                   "unpolarized": "unpolarized"}[polarization]
            ax.set_title(f"{material.name}: {pol} absorptance "
                         f"A(λ, θ) — Fresnel calculation")
            fig.tight_layout()
            return self._finalize(fig, filename), (lam_pk, ang_pk, A_pk)
