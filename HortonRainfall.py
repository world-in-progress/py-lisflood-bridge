import numpy as np
import math
from datetime import datetime, timedelta
from anuga import Domain

class HortonRainfall:
    """
    ANUGA-compatible rainfall + infiltration model
    """

    def __init__(
        self,
        times,            # seconds since start
        rain_rates,       # rainfall intensity (m/s)
        landuse,          # ndarray of landuse type per triangle
        start_time=0.0
    ):
        self.times = np.asarray(times)
        self.rain_rates = np.asarray(rain_rates)
        self.landuse = landuse
        self.start_time = start_time

        self.current_index = 0
        self.cumulative_rain_time = 0.0

    # ------------------------------------------------
    # Horton decay (same as your solver)
    # ------------------------------------------------
    @staticmethod
    def horton_decay(f0, fc, k, t):
        return fc + (f0 - fc) * math.exp(-k * t)

    def infiltration_rate(self, t):
        """
        Compute infiltration rate (m/s) for each triangle
        """
        infil = np.zeros_like(self.landuse, dtype=float)

        # time in hours
        t_hr = t / 3600.0

        # Horton params (inch/hr)
        params = {
            1: (0.8, 0.02, 10.0),   # building
            2: (0.8, 0.02, 10.0),   # road
            3: (3.0, 0.1, 2.0),     # farmland
            4: (0.0, 0.0, 0.0),     # pond
            5: (3.0, 0.1, 2.0),     # mountain
            6: (0.0, 0.0, 0.0),     # water
            7: (3.0, 0.1, 2.0),     # catch basin
        }

        for lu, (f0, fc, k) in params.items():
            mask = self.landuse == lu
            if f0 == 0.0:
                continue
            infil_in_hr = self.horton_decay(f0, fc, k, t_hr)
            infil[mask] = infil_in_hr * 0.0254 / 3600.0  # → m/s

        return infil

    # ------------------------------------------------
    # ANUGA interface
    # ------------------------------------------------
    def __call__(self, domain:Domain,t):
        """
        Return effective rainfall (m/s) at time t
        """
        t_rel = t - self.start_time

        # update rainfall index
        while (
            self.current_index + 1 < len(self.times)
            and t_rel > self.times[self.current_index + 1]
        ):
            self.current_index += 1

        rain = self.rain_rates[self.current_index]

        if rain > 0.0:
            self.cumulative_rain_time += domain.get_timestep()
            infil = self.infiltration_rate(self.cumulative_rain_time)
            return np.maximum(rain - infil, 0.0)

        return 0.0
