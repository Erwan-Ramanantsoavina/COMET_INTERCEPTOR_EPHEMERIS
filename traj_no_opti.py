import numpy as np
from godot.core import num, tempo, astro, events, ipfwrap
from godot.core import autodif as ad
from godot.core.autodif import bridge as br
import godot.core.util as c_util
from godot.model import interface, frames, common, prop
from godot import cosmos
from godot.cosmos import util
import numpy as np
import matplotlib.pyplot as plt
from ruamel import yaml
import time, os, copy
import pygmo as pg
import midas
from transfer_ephe import *
from totalDV import *

c_util.suppressLogger()

config_uni = util.load_yaml('Universe/universe.yml')
universe = cosmos.Universe( config_uni )


x0_halo_cr3bp = np.array([ 1.00737013  , 0.    ,     -0.00292083 , 0.     ,     0.01298799 , 0.        ])
T_halo_cr3bp = 3.0860603629916086

x1 = np.array([ 1.01102399e+00 , 8.48613753e-05 , 3.95486800e-03 , -1.63511915e-04 , -1.09131726e-02 , -1.28616305e-03])
x2 = np.array([ 1.00002001e+00 , -1.90647573e-04 , 7.04645004e-04 , 8.73703667e-02 , -2.38632662e-02 , -5.47251965e-02])
t1 = 4.157809328262431
t2 = 5.374938412849788
m= 0.49745782300811253

date_fly_by = tempo.Epoch('2025-10-22T04:48:00.0 TDB')

n_pt_traj = 3





# Computation

# Getting the right dates

tof_secondes = (t1 + t2) * T_SUN_EARTH/2/np.pi
date_start_from_halo = date_fly_by - tof_secondes

period_secondes = T_halo_cr3bp * T_SUN_EARTH/2/np.pi
date_start = date_start_from_halo - m * period_secondes

n_pt = 4
n_orb = 1

dv_halo = totalDV_halo(n_pt*n_orb , universe)
universe.evaluables.add('TotalDV_halo' , dv_halo)

dx = 100000
dv = 0.2
scales_vect = [ 2000000, 870000 , 600000, dv, dv, dv]
delta_vect = np.array([dx, dx, dx, dv , dv , dv])

config_traj , config_prob = config_halo(x0_halo_cr3bp , T_halo_cr3bp , date_start , n_pt , n_orb , scales_vect , delta_vect)



traj = cosmos.Trajectory(universe, config_traj)

prob = cosmos.Problem(universe, [traj], config_prob, useGradient=True)
# plot initial guess
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

x = universe.frames.vector6(
    'Earth',        # corps central
    'SC_center',    # point du spacecraft
    'SEROT',        # repère
    date_start_from_halo               # epoch
)


conf_traj2 , conf_prob , conf_prob2 = config_trajectory(x , date_start_from_halo , x1 , x2 , t1 , t2 , dv , n_pt_traj , False)

config_uni = util.load_yaml('Universe/universe.yml')
universe = cosmos.Universe( config_uni )
dv_traj = totalDV_traj(universe)
universe.evaluables.add('TotalDV_traj' , dv_traj)
traj2 = cosmos.Trajectory(universe, conf_traj2)

prob2 = cosmos.Problem(universe, [traj2], conf_prob2, useGradient=True)


traj2.compute(partials=False)

util.save_yaml(conf_traj2, "traj_example.yml")
util.save_yaml(conf_prob2, "problem.yml")
fig, ax = plt.subplots(figsize=(10, 6))  # Rectangle horizontal, haute résolution

ax.set_aspect('equal', adjustable='box')
ax.grid(True)
ax.set_xlabel(r'$X$ (km)')
ax.set_ylabel(r'$Y$ (km)')
ax.set_title("Intertial XY")

# Grille temporelle
E_grid = tempo.EpochRange(traj2.point('ctr0a'), traj2.point('end_pointa')).createGrid(3600)

XROT = []
for E_ in E_grid:
    x_rot = universe.frames.vector6('Earth', 'SC_center', 'SEROT', E_)
    XROT.append(np.concatenate(([E_.mjd()], x_rot)))

XROT = np.array(XROT)

# Temps des manoeuvres
t_m2 = traj2.point('ctr2a')

# Séparation des segments
seg1 = XROT[XROT[:, 0] <= t_m2.mjd()]   # avant 2e manoeuvre
seg2 = XROT[XROT[:, 0] >  t_m2.mjd()]   # après 2e manoeuvre

# Tracés
ax.plot(seg1[:, 1], seg1[:, 2], color='red', label='Manifold arc')
ax.plot(seg2[:, 1], seg2[:, 2], color='blue' , label='Lambert arc')

# Point final (intercept)
x_final, y_final = XROT[-1, 1], XROT[-1, 2]
ax.scatter(x_final, y_final, marker='x', s=50,  color='red', label='Intercept point', zorder=5)

# Terre
ax.scatter(0, 0, s=50, color='orange', label='Earth', zorder=5)


cp_names = [f'ctr{i}a' for i in range(2*n_pt_traj)]

for cp in cp_names:
    P = traj2.point(cp)

    x_cp = universe.frames.vector6('Earth', 'SC_center', 'SEROT', P)

    ax.scatter(
        x_cp[0], x_cp[1],
        s=25,
        color='green',
        marker='o',
        zorder=6
    )

#plt.tight_layout()

plt.show()
