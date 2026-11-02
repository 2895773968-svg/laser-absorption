# -*- coding: utf-8 -*-
"""
generate_spectrum.py —— README 用正入射光谱图
==============================================
对比 Cu/Al/Fe 三种金属在垂直入射下的吸收率-波长曲线（发表级风格）。
运行：python src/generate_spectrum.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from material_optics import MaterialOptics
from plotter import Plotter

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
DATASETS = [("Cu", "Cu_Johnson.yml"),
            ("Al", "Al_Rakic.yml"),
            ("Fe", "Fe_Johnson.yml")]


def main():
    plotter = Plotter(style="publication",
                      out_dir=os.path.join(os.path.dirname(__file__), "..", "output"))
    materials = [MaterialOptics.from_yaml(name, os.path.join(DATA_DIR, fname))
                 for name, fname in DATASETS]
    path = plotter.plot_spectrum(materials, polarization="unpolarized",
                                 filename="spectrum_unpolarized")
    print(f"已生成: {path}")


if __name__ == "__main__":
    main()
