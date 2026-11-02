# -*- coding: utf-8 -*-
"""让 pytest 能找到 src/ 下的模块（无需安装为包的轻量方案）。"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
