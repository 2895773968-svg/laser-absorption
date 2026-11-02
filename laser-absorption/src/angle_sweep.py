# -*- coding: utf-8 -*-
"""
angle_sweep.py —— AngleSweep 类
================================
固定波长，扫描入射角，输出 s/p/非偏振三种偏振态的吸收率曲线。

设计说明：
    AngleSweep 是 MaterialOptics 的"查询视图"——它不持有物理公式，
    全部计算委托给 MaterialOptics.absorptance()；
    结果封装为不可变的数据类 AngleSweepResult，供 Plotter 与测试直接消费。
"""

from dataclasses import dataclass, field
import numpy as np

from material_optics import MaterialOptics


@dataclass(frozen=True)
class AngleSweepResult:
    """一次角度扫描的完整结果（frozen=True 防止意外篡改物理数据）。"""
    material_name: str
    wavelength_nm: float
    angles_deg: np.ndarray          # 形状 (Nθ,)
    A_s: np.ndarray                 # s 偏振吸收率
    A_p: np.ndarray                 # p 偏振吸收率
    A_un: np.ndarray                # 非偏振（平均）吸收率
    reference: str = field(default="")

    @property
    def p_peak(self):
        """p 偏振吸收率最大值及其对应角度（金属'伪布儒斯特角'）。

        物理含义：p 偏振在掠入射区出现吸收峰，
        是金属表面反射系数相位变化导致的类布儒斯特效应。
        """
        i = int(np.argmax(self.A_p))
        return self.A_p[i], float(self.angles_deg[i])


class AngleSweep:
    """对某种金属在固定波长下做入射角扫描。"""

    def __init__(self, material: MaterialOptics, wavelength_nm: float):
        """
        参数
        ----
        material     : MaterialOptics 实例（数据源）
        wavelength_nm: 固定的激光波长，如 1064.0
        """
        self.material = material
        self.wavelength_nm = float(wavelength_nm)

    def run(self, angle_range=(0.0, 89.0), num=178) -> AngleSweepResult:
        """执行扫描。

        角度上限取 89 度而非 90 度：θ→90° 时 cosθ→0，
        菲涅尔系数趋于奇点（数值上仍有限，但物理上掠入射光斑发散，
        不属于有效加工区间），毕设中建议显式避开。
        """
        angles = np.linspace(angle_range[0], angle_range[1], num)
        m = self.material
        A_s = m.absorptance(self.wavelength_nm, angles, "s")
        A_p = m.absorptance(self.wavelength_nm, angles, "p")
        A_un = m.absorptance(self.wavelength_nm, angles, "unpolarized")
        return AngleSweepResult(m.name, self.wavelength_nm, angles,
                                A_s, A_p, A_un, m.reference)
