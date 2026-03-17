import numpy as np
from godot.core import num, tempo, astro, events, ipfwrap
from godot.core import autodif as ad
from godot.core.autodif import bridge as br
import godot.core.util as c_util
from godot.model import interface, frames, common, prop
from godot import cosmos
from godot.cosmos import util
import numpy as np
import matplotlib.pyplot as plt # Plotting like maltab
from ruamel import yaml
import time, os, copy
import pygmo as pg
import pygmo_plugins_nonfree as ppnf
import midas
from transfer_ephe import *
from totalDV import *


c_util.suppressLogger()

config_uni = util.load_yaml('Universe/universe.yml')
universe = cosmos.Universe( config_uni )

x0_cr3bp = [ 1.00737013  , 0.    ,     -0.00292083 , 0.     ,     0.01298799 , 0.        ]
T_cr3bp = 3.0860603629916086

date_start = tempo.Epoch('2029-01-01T00:00:00.0 TDB')
date_end = tempo.Epoch('2029-12-01T00:00:00.0 TDB')

n_pt = 4
n_orb = 3

dv_halo = totalDV_halo(n_pt*n_orb , universe)
universe.evaluables.add('TotalDV_halo' , dv_halo)
dx = 100000
dv = 0.2
scales_vect = [ 2000000, 870000 , 600000, dv, dv, dv]
delta_vect = np.array([dx, dx, dx, dv , dv , dv])

config_traj , config_prob = config_halo(x0_cr3bp , T_cr3bp , date_start , n_pt , n_orb , scales_vect , delta_vect)


traj = cosmos.Trajectory(universe, config_traj)

prob = cosmos.Problem(universe, [traj], config_prob, useGradient=True)
# plot initial guess
fig,ax = plt.subplots(1,3)
ax[0].grid()
ax[0].set_aspect('equal')
ax[0].set_xlabel('$X_{rot}$ (km)')
ax[0].set_ylabel('$Y_{rot}$ (km)')
ax[1].grid()
ax[1].set_aspect('equal')
ax[1].set_xlabel('$X_{rot}$ (km)')
ax[1].set_ylabel('$z_{rot}$ (km)')
ax[2].grid()
ax[2].set_xlabel('$Y_{rot}$ (km)')
ax[2].set_ylabel('$z_{rot}$ (km)')

traj.compute(partials=False)

# Prepare Pygmo problem
problem = pg.problem(prob)
tol_con = 1e-5
problem.c_tol = [tol_con] * problem.get_nc()


# Prepare population with initial guess
x0 = prob.get_x()
pop = pg.population(problem, 0)
pop.push_back(x0)

# Solve!
ip = pg.ipopt()
ip.set_numeric_option("tol",1e-3)
algo = pg.algorithm(ip)
algo.set_verbosity(1)
pop = algo.evolve(pop)

#################################################################################
# Solution
traj.compute(partials=False)
traj_up = traj.applyParameterChanges()
conf_update = util.deep_update(config_traj, traj_up)


dv_tot = universe.evaluables.get('SC_dv').eval(traj.point('end_point'))
print(f'Total dv: {dv_tot*1e3} m/s')
t_tot = traj.point('end_point') - traj.point('ctr0')
print(f'Total time: {t_tot/86400} days')


sol = traj.getTimelineSolution()
dv_tot = 0
for li in sol :
    for t in li :
        if "man" in t.name and (not t.name.endswith("end") and not t.name.endswith("start")):
            dvx = universe.evaluables.get(t.name+'_dv_x').eval(t.epoch)
            dvy = universe.evaluables.get(t.name+'_dv_y').eval(t.epoch)
            dvz = universe.evaluables.get(t.name+'_dv_z').eval(t.epoch)
            dv_tot += np.sqrt(dvx**2 + dvy**2 + dvz**2)

print(dv_tot)
E_grid = tempo.EpochRange(traj.point('ctr0'),traj.point('end_point')).createGrid(86400)
XROT=[]
dv_g = []
epog = []
for E_ in E_grid:
    x_rot = universe.frames.vector6('Earth','SC_center', 'SEROT', E_)
    XROT.append(np.concatenate(([E_.mjd()],x_rot)))
    epog.append(E_.mjd())
    dv_g.append(universe.evaluables.get('SC_dv').eval(E_))

ax[0].plot([x_[1] for x_ in XROT],[x_[2] for x_ in XROT] , color = 'blue' , linewidth = 0.3)
ax[1].plot([x_[1] for x_ in XROT],[x_[3] for x_ in XROT] , color = 'blue' , linewidth = 0.3)
#ax[2].plot([x_[2] for x_ in XROT],[x_[3] for x_ in XROT] , color = 'blue' , linewidth = 0.3)
ax[2].plot(epog , dv_g)

f = totalDV_halo(n_pt*n_orb , universe)
print(f.eval(traj.point('end_point')))
print(universe.evaluables.get('TotalDV_halo').eval(traj.point('end_point')))
plt.show()
