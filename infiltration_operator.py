import numpy as np
from anuga import Operator, Domain, Quantity

import numpy as np
import anuga

INCH_TO_M = 0.0254

def horton_decay(initial, final, k, t_hr):
    return final + (initial - final) * np.exp(-k * t_hr)


class HortonInfiltrationByTypeOperator(Operator):
    """
    Horton infiltration operator based on landuse type per triangle
    """

    def __init__(self, domain:Domain, tri_types):
        super().__init__(domain)

        self.tri_types = np.asarray(tri_types, dtype=np.int32)

        self.stage: Quantity = domain.quantities['stage']
        self.elev: Quantity = domain.quantities['elevation']

        # 累计“正在降雨”的时间（秒）
        self.cumulative_rain_time = 0.0
        self.last_time = domain.get_time()

        # Horton 参数表（inch/hr）
        self.params = {
            1: (0.8, 0.02, 10.0),
            2: (0.8, 0.02, 10.0),
            3: (3.0, 0.1, 2.0),
            5: (3.0, 0.1, 2.0),
            7: (3.0, 0.1, 2.0),
            # 4,6 默认 0
        }
        self.frs: list[float] = [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        ]

        # Initialize frs
        for land_type, (f0, fc, k) in self.params.items():
            fr_index = land_type - 1
            self.frs[fr_index] = horton_decay(f0, fc, k, 0.0) * INCH_TO_M / 3600.0

    def __call__(self, domain:Domain=None):
        if domain is None:
            domain = self.domain

        t = domain.get_time()
        dt = t - self.last_time
        if dt <= 0.0:
            return

        self.last_time = t

        # -------------------------
        # 判断是否在下雨
        # -------------------------
        rain_rate = 0.0
        for forcing in domain.forcing_terms:
            if isinstance(forcing, anuga.Rainfall):
                rain_rate = forcing.rate(t) if callable(forcing.rate) else forcing.rate
                break
        
        # Update frs
        if rain_rate > 0.0:
            self.cumulative_rain_time += dt
            t_hr = self.cumulative_rain_time / 3600.0

            for land_type, (f0, fc, k) in self.params.items():
                fr_index = land_type - 1
                self.frs[fr_index] = horton_decay(f0, fc, k, t_hr) * INCH_TO_M / 3600.0

        # -------------------------
        # 逐 type 处理
        # -------------------------
        elev = self.elev.centroid_values
        stage = self.stage.centroid_values
        for land_type, (f0, fc, k) in self.params.items():
            mask = (self.tri_types == land_type)
            if not np.any(mask):
                continue
            
            # 扣减水深（不能低于地面）
            fr_index = land_type - 1
            stage[mask] -= self.frs[fr_index] * dt
            stage[mask] = np.maximum(stage[mask], elev[mask])

        # water body / fish pond（4,6）自动为 0，不处理


# class TestOperator(Operator):
#     def __init__(self, domain: Domain, f0, fc, k, description = None, label = None, logging: bool = False):
#         Operator.__init__(self, domain, description, label, logging)

#         self.f0 = f0 / 1000.0 / 3600.0
#         self.fc = fc / 1000.0 / 3600.0
#         self.k = k / 3600.0

#         self.domain: Domain

#         self.start_time = 0.0
    
#     def __call__(self):
#         print('hahaha')
#         timestep = self.domain.get_timestep()
#         current_time = self.domain.get_time()

#         t = current_time
#         rate = self.fc + (self.f0 - self.fc) * np.exp(-self.k * t)
#         infiltration_depth = rate * timestep

#         # 3. 获取当前水深 (stage - elevation)
#         quantities = self.domain.quantities
#         stage: np.array = quantities['stage'].get_values()
#         elevation = quantities['elevation'].get_values()
#         depth = stage - elevation
        
#         # 4. 执行下渗：不能吸干水 (下渗量不能超过当前水深)
#         # 创建一个更新数组
#         actual_infiltration = np.minimum(depth, infiltration_depth)
        
#         # 只在有水的地方下渗 (depth > 0)
#         actual_infiltration = np.maximum(actual_infiltration, 0.0)
        
#         # 5. 更新 Stage (水位)
#         new_stage = stage - actual_infiltration
#         quantities['stage'].set_values(new_stage)


# if __name__ == '__main__':
#     a = TestOperator()
