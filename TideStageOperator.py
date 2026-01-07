import numpy as np
from anuga import Operator, Domain,Quantity

import numpy as np
from anuga import Operator

class TideStageOperator(Operator):
    """
    强制指定一批单元的 stage = elevation + tide(t)
    兼容 ANUGA 两种调用方式： operator() 或 operator(domain)
    """

    def __init__(self, domain:Domain, cell_indices, tide_func):
        super().__init__(domain)
        # 接受 list/tuple/ndarray
        self.cell_indices = np.asarray(cell_indices, dtype=int)
        self.tide_func = tide_func

        # 直接持有 quantity 引用（方便更新）
        # 声明类型为 Quantity
        self.stage: Quantity = domain.quantities['stage']
        self.elev: Quantity = domain.quantities['elevation']

    def __call__(self, domain:Domain=None):
        """
        兼容签名：ANUGA 可能会调用 operator() 或 operator(domain)
        如果 domain 参数为 None，就使用 self.domain（Operator 基类中保存的 domain）
        """
        if domain is None:
            domain = self.domain

        t = domain.get_time()
        tide = self.tide_func(t)

        # 直接在 centroid_values 上做批量赋值（保证索引为整数 ndarray）
        idx = self.cell_indices
        # 可选：防护索引越界（仅在调试需要）
        # idx = idx[idx >= 0]
        # idx = idx[idx < len(self.stage.centroid_values)]

        # 强制覆盖 stage
        self.stage.centroid_values[idx] = self.elev.centroid_values[idx] + tide

        # 可选调试输出（上线时可以注释掉）
        # print(f"[TideOp] t={t:.2f}, tide={tide:.4f}, n_cells={len(idx)}")
