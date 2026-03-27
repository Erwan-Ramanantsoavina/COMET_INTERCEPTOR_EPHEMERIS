from godot.core import tempo
from godot import cosmos
from godot.cosmos import util
import pygmo as pg
import godot.core.util as c_util
import numpy as np
import matplotlib.pyplot as plt

c_util.suppressLogger()

config_uni = util.load_yaml('Universe/universe.yml')# Load Universe
config_traj = util.load_yaml('traj_example.yml')    # Load Trajectory
config_prob = util.load_yaml('problem.yml')         # Load Problem

universe = cosmos.Universe(config_uni)              # Create Universe object
traj = cosmos.Trajectory(universe , config_traj)    # Create Trajectory object
udp = cosmos.Problem(universe , [traj] , config_prob , True)    # Create Problem object (Godot)

problem = pg.problem(udp)                           # Create problem object (Pygmo)
c_tol = 1e-6
problem.c_tol = [c_tol] * problem.get_nc()          # Set a tolerance for constraints

pop = pg.population(problem , 0)                    # initiate a population of size 0
x0 = udp.get_x()                                    # Get the initial decision vector
pop.push_back(x0)                                   # Put it in the population

algo = pg.algorithm(pg.ipopt())                     # Use the algorithm ipopt for example
pop = algo.evolve(pop)                              # Solve the problem

traj.applyParameterChanges()                        # Apply the Changes from the optimization
# Plotting

epoch_start = traj.point('ctr0a')                   # Starting date of the trajectory
epoch_end = traj.point('end_pointa')                # ending date of the trajectory
step = 1 * 24 * 60 * 60                             # Step for the grid (1 day)

epoch_grid = tempo.EpochRange(epoch_start , epoch_end).createGrid(step) # create a grid of epochs

x_list = []
dv_list = []
epoch_mjd_list = []
for e in epoch_grid :
    x = universe.frames.vector6('Earth' , 'SC_center' , 'SEROT' , e)    # State vector of SC_center from the Earth, considering the axes SEROT at the date e
    dv = universe.evaluables.get('SC_dv').eval(e)                       # dv at the date e
    x_list.append(x)
    dv_list.append(dv*1000)
    epoch_mjd_list.append(e.mjd())                                      # to have a list od epoch in mjd2000 format

x_list = np.array(x_list)

dv_final = universe.evaluables.get('SC_dv').eval(epoch_end)             # Get the Delta-v total final
print(f"Final DV = {dv_final} km/s")
plt.figure()
plt.grid()
plt.axis('equal')
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.plot(x_list[:,0] , x_list[:,1])

plt.figure()
plt.grid()
plt.xlabel('Epoch (mjd)')
plt.ylabel('Delta-v (m/s)')
plt.plot(epoch_mjd_list , dv_list)
plt.show()

    

