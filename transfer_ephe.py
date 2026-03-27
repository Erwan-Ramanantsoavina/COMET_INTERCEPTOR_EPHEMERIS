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
from aux_fun import ctr , man , match

D_SUN_EARTH = 149597870
GM_SUN = 132712440018	
GM_EARTH = 398600.4418
mu_se = GM_EARTH / (GM_EARTH + GM_SUN)
T_SUN_EARTH = 2*np.pi*np.sqrt(D_SUN_EARTH**3/(GM_EARTH+GM_SUN))

def config_halo(x0_cr3bp , T_CR3BP , date_start , n_pt , n_orb , scales_vect , delta_vect) :
    
    CRTBP = midas.astro.crtbp.AdimensionalCRTBP(mu_se)
    IntegratorCRTBP = midas.astro.crtbp.AdimensionalIntegratorCRTBP(mu_se)
    config_uni = util.load_yaml('Universe/universe.yml')
    universe = cosmos.Universe( config_uni )
    period = T_CR3BP * T_SUN_EARTH/(2*np.pi*86400)
    dt = period / n_pt

    mjd_vect_ref = date_start.mjd() + np.arange(0, n_orb*period + dt, dt)

    x0_vect_ref = []

    for m_ in mjd_vect_ref:

        # Get CR3BP state
        dt_cr3bp = (m_ - date_start.mjd()) % period
        dt_cr3bp = dt_cr3bp*86400*2*np.pi/T_SUN_EARTH
        [t_cr3bp, x_cr3bp] = IntegratorCRTBP.integrate(x0_cr3bp, dt_cr3bp)
        x0_rot = np.copy(x_cr3bp[-1])
        
        lu = universe.frames.distance('Sun','Earth', tempo.Epoch(str(m_) + ' TDB'))
        tu = np.sqrt(lu**3/(GM_SUN+GM_EARTH))
        x0_rot[0] -= 1 - mu_se
        x0_rot[:3] *= lu
        x0_rot[3:] *= lu/tu
        x0_vect_ref.append(x0_rot)
        

    #################################################################################
    # Prepare timeline
    config_traj = {'settings': {'relTol': 1e-12, 'steps': 10000000},
               'setup': [
                   {'name': 'SC',
                    'type': 'group', 
                    'spacecraft' : 'SC',
                    'input': [
                    {'name': 'center',
                        'type': 'point'},
                    {'name': 'dv',
                        'type': 'scalar',
                        'unit': 'm/s'}
                        ]}]
                    }
    TT = []
    i_ = 0
    x_point_ref = []
    # Add points
    for m_, x_ in zip(mjd_vect_ref, x0_vect_ref):
        ct = ctr('ctr'+str(i_) , tempo.Epoch(str(m_) + ' TDB') , x_ )
        TT.append(ct)
        
        if i_ < n_orb*n_pt:

            mn = man('man'+str(i_) , 'ctr'+str(i_) , 0 , np.array([0.0001 , 0.0001 , 0.0001]))
            mtch = match('match'+str(i_) , 'ctr'+str(i_) , dt*86400/2 , 'ctr'+str(i_+1) , dt*86400/2)
            TT.append(mn)
            TT.append(mtch)

        i_ += 1
        x_point_ref.append(x_)

    # add final point
    end_point = {'type': 'point',
             'input': 'SC',
             'name': 'end_point',
             'point': {'reference': 'ctr'+str(i_ - 1),
                       'dt': '1.0 s'}}
    TT.append(end_point)

    # Assemble trajectory
    config_traj['timeline'] = TT

    #################################################################################
    # Assemble problem
    config_prob = {
    'objective': {'type': 'minimise',
                  'point': 'end_point',
                  'value': 'TotalDV_halo',
                  'scale': 1e6},
    'parameters': {'free': []}}
    FF = []
    SS = {}
    BB = {}
    cart_vect = ['pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z']

    units_vect = [' km', ' km', ' km', ' km/s', ' km/s', ' km/s']

    
    for i_ in range(n_pt*n_orb):
        if not (i_==0): 
            FF.append('ctr' + str(i_)+'_dt')
            FF.append('ctr' + str(i_)+'_SC_dv')
            SS['ctr' + str(i_)+'_dt'] = '1 day'
            SS['ctr' + str(i_)+'_SC_dv'] = '1 mm/s'
            BB['ctr' + str(i_)+'_dt'] = ['-1 day', '1 day']
            BB['ctr' + str(i_)+'_SC_dv'] = ['0 m/s', '1000 m/s']
        for c_, s_, u_, b_, d_ in zip(cart_vect, scales_vect, units_vect, x_point_ref[i_],delta_vect):
            FF.append('ctr'+str(i_)+'_SC_center_'+c_)
            SS['ctr'+str(i_)+'_SC_center_'+c_] = str(s_)+u_
            BB['ctr'+str(i_)+'_SC_center_'+c_] = [str(b_-d_)+u_,str(b_+d_)+u_]

        # Add manoeuvre
        FF.append('man'+str(i_)+'_dv_x')
        FF.append('man'+str(i_)+'_dv_y')
        FF.append('man'+str(i_)+'_dv_z')

        SS['man'+str(i_)+'_dv_x'] = '1 mm/s'
        SS['man'+str(i_)+'_dv_y'] = '1 mm/s'
        SS['man'+str(i_)+'_dv_z'] = '1 mm/s'
        BB['man'+str(i_)+'_dv_x'] = ['-1 m/s', '1 m/s']
        BB['man'+str(i_)+'_dv_y'] = ['-1 m/s', '1 m/s']
        BB['man'+str(i_)+'_dv_z'] = ['-1 m/s', '1 m/s']
        # Add match point
        FF.append('match'+str(i_)+'_right_dt')
        SS['match'+str(i_)+'_right_dt'] = '1 day'
        BB['match'+str(i_)+'_right_dt'] = ['1 s', str(period/n_pt+1)+' day']
    # Add last point
    i_ = n_pt*n_orb 
    FF.append('ctr' + str(i_)+'_dt')
    FF.append('ctr' + str(i_)+'_SC_dv')
    SS['ctr' + str(i_)+'_dt'] = '1 day'
    SS['ctr' + str(i_)+'_SC_dv'] = '1 mm/s'
    BB['ctr' + str(i_)+'_dt'] = ['-1 day', '1 day']
    BB['ctr' + str(i_)+'_SC_dv'] = ['0 m/s', '100 m/s']
    for c_, s_, u_, b_, d_ in zip(cart_vect, scales_vect, units_vect, x_point_ref[0],delta_vect):
        FF.append('ctr'+str(i_)+'_SC_center_'+c_)
        SS['ctr'+str(i_)+'_SC_center_'+c_] = str(s_)+u_
        BB['ctr'+str(i_)+'_SC_center_'+c_] = [str(b_-d_)+u_,str(b_+d_)+u_]

    
    

    # Setup bounds on left and right points
    config_prob['parameters']['free'] = FF
    config_prob['parameters']['scales'] = SS
    config_prob['parameters']['bounds'] = BB
    
    return config_traj , config_prob

def config_trajectory(x0_halo_ephe , date_start , x1_cr3bp , x2_cr3bp , t1_cr3bp , t2_cr3bp , dv_cr3bp , n_pt , with_man = False) :
        
    # create trajectory
    config_traj = {'settings': {'relTol': 1e-12, 'steps': 10000000},
               'setup': [
                   {'name': 'SC',
                    'type': 'group', 
                    'spacecraft' : 'SC',
                    'input': [
                    {'name': 'center',
                        'type': 'point'},
                    {'name': 'dv',
                        'type': 'scalar',
                        'unit': 'm/s'}
                        ]}]
                    }
    
    # TWO problems :
    # -Assembling the arcs
    # -Assembling the arcs while minimizing th DV.
    config_prob = {'parameters': {'free': []}}
    
    config_prob2 = {
        'objective' :{'type' : 'minimise' ,
                      'point' : 'end_pointa' ,
                      'value' : 'TotalDV_traj' ,
                      'scale' : 1e3} ,

        'parameters': {'free': []}}


    list_ctr = []
    FF = []
    FF2 = []

    # Create SUN-EARTH CR3BP
    CRTBP = midas.astro.crtbp.AdimensionalCRTBP(mu_se)
    IntegratorCRTBP = midas.astro.crtbp.AdimensionalIntegratorCRTBP(mu_se)

    # Load Universe
    uni_config = cosmos.util.load_yaml('Universe/universe.yml')
    uni = cosmos.Universe(uni_config)

    # List for the 'timeline' parameter
    TT = [] 

    # Add the Starting Point on the Halo orbit.
    ctr_departure = ctr('ctr0a' , date_start , x0_halo_ephe)
    TT.append(ctr_departure)
    list_ctr.append(x0_halo_ephe)

    # Add a maneuver at the start : injection on the Unstable direction.
    lu = uni.frames.distance('Sun','Earth', date_start)
    tu = np.sqrt(lu**3/(GM_SUN+GM_EARTH))
    dv1_vect_kms = x1_cr3bp[3:]*lu/tu - x0_halo_ephe[3:]
    man_injection = man('man_injection' , 'ctr0a' , 0 , dv1_vect_kms)
    TT.append(man_injection)

    # Add control and match points with a dt fixed
    t1_s = t1_cr3bp * np.sqrt(D_SUN_EARTH**3/(GM_SUN+GM_EARTH))
    dt1_s = t1_s / n_pt
    
    match0 = match('match0a' , 'ctr0a' , dt1_s/2 , 'ctr1a' , dt1_s/2)
    TT.append(match0)
    
    for i in range(1 , n_pt+1) :

        # Create the ctr point
        [t_cr3bp , x_cr3bp] = IntegratorCRTBP.integrate(x1_cr3bp , i*t1_cr3bp/n_pt)
        x_loc = x_cr3bp[-1]
        lu = uni.frames.distance('Sun','Earth', date_start + i*dt1_s)
        tu = np.sqrt(lu**3/(GM_SUN+GM_EARTH))
        x_loc[0] -= 1-mu_se
        x_loc[:3] *= lu
        x_loc[3:] *= lu/tu
        list_ctr.append(x_loc)
        ctr_intermediaire = ctr('ctr'+str(i)+'a' , date_start + i*dt1_s , x_loc)
        TT.append(ctr_intermediaire)

        if i < n_pt :
            if with_man :
                mani = man('man'+str(i)+'a' , 'ctr'+str(i)+'a' , 0 , np.array([0.00001 , 0.00001 , 0.00001]))
                TT.append(mani)
            matchi = match('match'+str(i)+'a' , 'ctr'+str(i)+'a' , dt1_s/2 , 'ctr'+str(i+1)+'a' , dt1_s/2)
            TT.append(matchi)
            


    # Add the dsm maneuver
    d = uni.frames.distance('Sun','Earth', date_start + t1_s)
    t = np.sqrt(lu**3/(GM_SUN+GM_EARTH))
    dsm_vect = x2_cr3bp[3:]*d/t - x_loc[3:]
    man_dsm = man('man_dsm' , 'ctr'+str(n_pt)+'a' , 0 , dsm_vect)
    TT.append(man_dsm)

    # Then again add the matching and control points
    t2_s = t2_cr3bp * np.sqrt(D_SUN_EARTH**3/(GM_SUN+GM_EARTH))
    dt2_s = t2_s/n_pt

    match2 = match('match'+str(n_pt)+'a' , 'ctr'+str(n_pt)+'a' , dt2_s/2 , 'ctr'+str(n_pt+1)+'a' , dt2_s/2)
    TT.append(match2)

    for i in range(1 , n_pt + 1) :
        [t_cr3bp , x_cr3bp] = IntegratorCRTBP.integrate(x2_cr3bp , i*t2_cr3bp/n_pt)
        x_loc = x_cr3bp[-1]
        lu = uni.frames.distance('Sun','Earth', date_start + t1_s + i*dt2_s)
        tu = np.sqrt(lu**3/(GM_SUN+GM_EARTH))
        x_loc[0] -= 1-mu_se
        x_loc[:3] *= lu
        x_loc[3:] *= lu/tu
        list_ctr.append(x_loc)
        ctr_intermediaire = ctr('ctr'+str(n_pt + i)+'a' , date_start + t1_s + i*dt2_s , x_loc)
        TT.append(ctr_intermediaire)

        if i < n_pt :
            if with_man :
                mani = man('man'+str(n_pt + i)+'a' , 'ctr'+str(n_pt + i)+'a' , 0 , np.array([0.00001 , 0.00001 , 0.00001]))
                TT.append(mani)
            match3 = match('match'+str(n_pt + i)+'a' , 'ctr'+str(n_pt+i)+'a' , dt2_s/2 , 'ctr'+str(n_pt+1+i)+'a' , dt2_s/2)
            TT.append(match3)
            



    end_point = {'type': 'point',
             'input': 'SC',
             'name': 'end_pointa',
             'point': {'reference': 'ctr'+str(2*n_pt)+'a',
                       'dt': '0 s'}}
    TT.append(end_point)
    
    config_traj['timeline'] = TT


    # Les variables à bouger pour résoudre le problème

    # Le vecteur dv1
    # FF.append('ctr0a_dt')
    
    FF.append('man_injection_dv_x')
    FF.append('man_injection_dv_y')
    FF.append('man_injection_dv_z')

    FF2.append('ctr0a_dt')
    #FF2.append('man_injection_dt')
    FF2.append('man_injection_dv_x')
    FF2.append('man_injection_dv_y')
    FF2.append('man_injection_dv_z')

    # Les positions des points le long de la trajectoire
    for i in range(1 , 2*n_pt) :
        FF.append('ctr'+str(i)+'a_dt')
        FF.append('ctr'+str(i)+'a_SC_dv')
        FF.append('ctr'+str(i)+'a_SC_center_pos_x')
        FF.append('ctr'+str(i)+'a_SC_center_pos_y')
        FF.append('ctr'+str(i)+'a_SC_center_pos_z')
        FF.append('ctr'+str(i)+'a_SC_center_vel_x')
        FF.append('ctr'+str(i)+'a_SC_center_vel_y')
        FF.append('ctr'+str(i)+'a_SC_center_vel_z') 

        FF2.append('ctr'+str(i)+'a_dt')
        FF2.append('ctr'+str(i)+'a_SC_dv')
        FF2.append('ctr'+str(i)+'a_SC_center_pos_x')
        FF2.append('ctr'+str(i)+'a_SC_center_pos_y')
        FF2.append('ctr'+str(i)+'a_SC_center_pos_z')
        FF2.append('ctr'+str(i)+'a_SC_center_vel_x')
        FF2.append('ctr'+str(i)+'a_SC_center_vel_y')
        FF2.append('ctr'+str(i)+'a_SC_center_vel_z')   

    # Le vecteur dv2
    FF.append('man_dsm_dv_x')
    FF.append('man_dsm_dv_y')
    FF.append('man_dsm_dv_z')

    #FF2.append('man_dsm_dt')
    FF2.append('man_dsm_dv_x')
    FF2.append('man_dsm_dv_y')
    FF2.append('man_dsm_dv_z')
    
    if with_man :
        for i in range(1,2*n_pt):
            if i != n_pt :
                FF.append('man'+str(i)+'a_dv_x')
                FF.append('man'+str(i)+'a_dv_y')
                FF.append('man'+str(i)+'a_dv_z')

                FF2.append('man'+str(i)+'a_dv_x')
                FF2.append('man'+str(i)+'a_dv_y')
                FF2.append('man'+str(i)+'a_dv_z')
  
    FF.append('ctr'+str(2*n_pt)+'a_SC_dv')       
    FF.append('ctr'+str(2*n_pt)+'a_SC_center_vel_x')
    FF.append('ctr'+str(2*n_pt)+'a_SC_center_vel_y')
    FF.append('ctr'+str(2*n_pt)+'a_SC_center_vel_z')

    # FF.append('man_injection_dt')
    # FF.append('man_dsm_dt')
    # FF2.append('man_injection_dt')
    # FF2.append('man_dsm_dt')

    FF2.append('ctr'+str(2*n_pt)+'a_SC_dv')       
    FF2.append('ctr'+str(2*n_pt)+'a_SC_center_vel_x')
    FF2.append('ctr'+str(2*n_pt)+'a_SC_center_vel_y')
    FF2.append('ctr'+str(2*n_pt)+'a_SC_center_vel_z')

    config_prob['parameters']['free'] = FF
    config_prob2['parameters']['free'] = FF2

    return config_traj , config_prob , config_prob2
