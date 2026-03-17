import godot.model.common as common
from godot.cosmos import util
import godot.core.util as c_util
from godot.core import autodif as ad
import numpy as np

c_util.suppressLogger()

class totalDV_traj( common.ScalarTimeEvaluable ) :

    def __init__(self, universe):
        super().__init__()
        self.universe = universe

    def eval(self, e ):
        

        dv_tot = 0.0

        dvx_inj = self.universe.evaluables.get(f'man_injection_dv_x').eval(e)
        dvy_inj = self.universe.evaluables.get(f'man_injection_dv_y').eval(e)
        dvz_inj = self.universe.evaluables.get(f'man_injection_dv_z').eval(e)
        dv_tot += ad.sqrt(dvx_inj*dvx_inj + dvy_inj*dvy_inj + dvz_inj*dvz_inj)

        dvx_dsm = self.universe.evaluables.get(f'man_dsm_dv_x').eval(e)
        dvy_dsm = self.universe.evaluables.get(f'man_dsm_dv_y').eval(e)
        dvz_dsm = self.universe.evaluables.get(f'man_dsm_dv_z').eval(e)
        dv_tot += ad.sqrt(dvx_dsm*dvx_dsm + dvy_dsm*dvy_dsm + dvz_dsm*dvz_dsm)
        return dv_tot


class totalDV_halo( common.ScalarTimeEvaluable ) :

    def __init__(self, n_man , universe):
        super().__init__()
        self.n_man = n_man
        self.universe = universe

    def eval(self, e ):
        

        dv_tot = 0.0
        for i__ in range(self.n_man) :
            dvx = self.universe.evaluables.get(f'man{i__}_dv_x').eval(e)
            dvy = self.universe.evaluables.get(f'man{i__}_dv_y').eval(e)
            dvz = self.universe.evaluables.get(f'man{i__}_dv_z').eval(e)
            dv_tot += ad.sqrt(dvx*dvx + dvy*dvy + dvz*dvz)

        return dv_tot

