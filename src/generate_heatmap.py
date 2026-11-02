# -*- coding: utf-8 -*-
"""
generate_heatmap.py —— README 主视觉图生成脚本
===============================================
功能：对 Cu/Al/Fe 各生成一张 "波长(400-1100nm) × 入射角(0-89°)"
的 p 偏振吸收率二维热力图，自动标注全局最大值点（伪布儒斯特吸收）。

运行：python src/generate_heatmap.py
输出：output/heatmap_<金属>_ppol.png / .pdf（300 DPI）
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # 支持直接运行

from material_optics import MaterialOptics
from plotter import Plotter

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

# (金属名, 数据文件) —— 与 data/ 目录中的 refractiveindex.info 文件对应
DATASETS = [("Cu", "Cu_Johnson.yml"),
            ("Al", "Al_Rakic.yml"),
            ("Fe", "Fe_Johnson.yml")]


def main():
    plotter = Plotter(style="publication",
                      out_dir=os.path.join(os.path.dirname(__file__), "..", "output"))

    print("=" * 66)
    print("p 偏振吸收率全局最大值（伪布儒斯特吸收的波长-角度位置）")
    print("-" * 66)
    for name, fname in DATASETS:
        mat = MaterialOptics.from_yaml(name, os.path.join(DATA_DIR, fname))
        path, (lam_pk, ang_pk, A_pk) = plotter.plot_heatmap(
            mat, polarization="p",
            filename=f"heatmap_{name}_ppol")
        print(f"{name:>3}: A_max = {A_pk:5.1f}%  @  λ = {lam_pk:6.1f} nm, "
              f"θ = {ang_pk:5.1f}°   ->  {os.path.basename(path)}")
    print("=" * 66)
    print("物理注解：短波长端（<600nm，铜的 d 带带间跃迁区）本底吸收率高，")
    print("与掠入射伪布儒斯特峰叠加，使最大值出现在左上角区域；")
    print("红外端（1064nm）本底仅 ~3%，但角度调制仍能把 A_p 拉高一个量级——")
    print("这正是'掠入射 + 短波长'双增益策略的物理基础。")


if __name__ == "__main__":
    main()
