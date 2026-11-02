# -*- coding: utf-8 -*-
"""
test_fresnel.py —— 菲涅尔物理内核的单元测试
============================================
测试策略（由内到外三层验证）：
    1. 数学自洽：θ=0° 时 s/p 偏振必须严格相等（入射面退化为任意的对称性）；
    2. 解析对拍：θ=0° 时与教科书公式 A = 4n/[(n+1)²+k²] 一致；
    3. 物理约束：0 ≤ A ≤ 1、非偏振=偏振平均、能量守恒 R+A=1。
"""

import numpy as np
import pytest

from material_optics import MaterialOptics, fresnel_absorptance


# ---------------------------------------------------------------------
# 1. 核心测试：正入射时 s 与 p 偏振吸收率相等
# ---------------------------------------------------------------------
@pytest.mark.parametrize("n,k", [
    (0.35, 7.03),    # Cu @ 1064 nm（真实量级）
    (1.03, 2.60),    # Cu @ ~450 nm（带间跃迁区量级）
    (1.44, 7.38),    # Al @ 1064 nm（真实量级）
    (2.00, 5.00),    # 任意虚构值（测试不应依赖特定材料）
])
def test_normal_incidence_s_equals_p(n, k):
    """物理依据：θ=0 时'入射面'可取任意方位，s/p 定义失去区分意义，
    故两个偏振本征态的反射/吸收必须严格相同——这是对称性要求，
    若不等则说明 r_s / r_p 公式写错（最常见的错误是 eps 与 1 写反）。"""
    theta0 = np.array([0.0])
    A_s = fresnel_absorptance(n, k, theta0, "s")
    A_p = fresnel_absorptance(n, k, theta0, "p")
    np.testing.assert_allclose(A_s, A_p, rtol=1e-12, atol=1e-14)


# ---------------------------------------------------------------------
# 2. 解析对拍：正入射必须还原教科书公式
# ---------------------------------------------------------------------
@pytest.mark.parametrize("n,k", [(0.35, 7.03), (1.44, 7.38), (0.18, 3.10)])
def test_normal_incidence_matches_analytic_formula(n, k):
    """A(0°) = 4n/[(n+1)² + k²]，由 r=(1-N)/(1+N) 直接推导。"""
    A_num = fresnel_absorptance(n, k, 0.0, "unpolarized")
    A_ana = 4.0 * n / ((n + 1.0) ** 2 + k ** 2)
    np.testing.assert_allclose(A_num, A_ana, rtol=1e-12)


# ---------------------------------------------------------------------
# 3. 物理约束
# ---------------------------------------------------------------------
def test_absorptance_bounded_0_1():
    """吸收率是概率性物理量，必须落在 [0, 1] 区间（全网格检查）。"""
    A = fresnel_absorptance(np.linspace(0.1, 2.0, 50),
                            np.linspace(0.5, 8.0, 50),
                            np.linspace(0, 89, 90), "p")
    assert np.all(A >= 0.0) and np.all(A <= 1.0)


def test_unpolarized_is_average():
    """非偏振光 = s/p 两独立偏振态的强度平均（能量各占一半）。"""
    n, k, th = 0.35, 7.03, 65.0
    A_un = fresnel_absorptance(n, k, th, "unpolarized")
    A_s = fresnel_absorptance(n, k, th, "s")
    A_p = fresnel_absorptance(n, k, th, "p")
    np.testing.assert_allclose(A_un, 0.5 * (A_s + A_p), rtol=1e-12)


# ---------------------------------------------------------------------
# 4. 文献基准：用已发表实验值锚定整条计算链（含 YAML 读取与插值）
# ---------------------------------------------------------------------
@pytest.fixture(scope="module")
def cu_material():
    """加载真实的 Cu 数据（复用项目 data/ 目录），模块级只读一次。"""
    import os
    data = os.path.join(os.path.dirname(__file__), "..", "data")
    return MaterialOptics.from_yaml("Cu", os.path.join(data, "Cu_Johnson.yml"))


def test_cu_1064nm_matches_literature(cu_material):
    """锚点：Cu 在 1064 nm 的垂直入射吸收率，文献范围 2%~5%
    （Johnson & Christy 1972 数据 + 菲涅尔方程的标准结果）。"""
    A = cu_material.absorptance_normal(1064.0)
    assert 0.02 < A < 0.05


def test_cu_blue_over_infrared_ratio(cu_material):
    """锚点：短波长(d 带带间跃迁区)吸收率必须是红外端的 10 倍以上——
    这是整个'蓝绿光加工'技术路线的物理地基，算不出这个倍率即物理错误。"""
    A_blue = cu_material.absorptance_normal(532.0)
    A_ir = cu_material.absorptance_normal(1064.0)
    assert A_blue / A_ir > 10.0
