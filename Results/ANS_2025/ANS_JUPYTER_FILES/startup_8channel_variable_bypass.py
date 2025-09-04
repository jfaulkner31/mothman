"""
This file is meant to model the MSRE startup transient
but with a modified bypass flow rate.

Runs through command line and saves results based on user input: startup_bypass_modified_flowrate_<PERCENT>_percent.pkl


Run this file using the following:
python startup_8channel_variable_bypass.py <PERCENT>
where <PERCENT> is the percent of flow that goes to the bypass flow channel.
"""

import openmc
import numpy as np
from Meshing.Meshing import *
from Fields.Fields import *
from Kernels.Kernels import *
from Solvers.Solvers import *
from Subchannel.FluidRelation import FluidRelation
from Subchannel.Channel import Channel
from Subchannel.Channel import ChannelInterface
from Aux.CSVObjects import *
from Aux.ReactorPhysicsObjects import *

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.patches import Rectangle

import pandas as pd

from Data.msre_data import msre_data_dict

import sys

"""
Setup
"""

# GET ARGS FROM INPUT

PERCENT = float(sys.argv[1])
print(f"Script name: {sys.argv[0]}")
print("BYPASS FLOW RATIO IS (%):", PERCENT)

# OUTPUT FILENAME
OUTPUT_FILENAME = "Results/ANS_2025/startup_bypass_modified_flowrate_"+str(PERCENT)+"_percent"+".pkl"

# FLUID
fluid = msre_data_dict['fluid']

### BC DATA ###
# mdot_csv = load_csv('Data/msre_flow_coastdown.csv')
mdot_csv = load_csv('Data/mass_flow_startup.csv')
MDOT_MAX = msre_data_dict['mdot_max']
START_TIME = -100
mdot_bc = MDOT_MAX / 100 * csv_interpolator(csv_df=mdot_csv, x_value=START_TIME, x_label='time', y_label='mdot')

T_bc = msre_data_dict['temp_bc']
pressure_bc = msre_data_dict['pressure_bc']

# SET CHANNELS UP FROM FULL CORE SOLUTION (NEED TO RUN COASTDOWN SIMULATION FIRST)
ch1 = Channel.import_from_pkl(filename='Data/model_8_channels/channel_1_setup.pkl')
ch2 = Channel.import_from_pkl(filename='Data/model_8_channels/channel_2_setup.pkl')
ch3 = Channel.import_from_pkl(filename='Data/model_8_channels/channel_3_setup.pkl')
ch4 = Channel.import_from_pkl(filename='Data/model_8_channels/channel_4_setup.pkl')
ch5 = Channel.import_from_pkl(filename='Data/model_8_channels/channel_5_setup.pkl')
ch6 = Channel.import_from_pkl(filename='Data/model_8_channels/channel_6_setup.pkl')
ch7 = Channel.import_from_pkl(filename='Data/model_8_channels/channel_7_setup.pkl')
ch8 = Channel.import_from_pkl(filename='Data/model_8_channels/channel_8_setup.pkl')

upper_plenum_channel = Channel.import_from_pkl(filename='Data/model_8_channels/upper_plenum.pkl')
lower_plenum_channel = Channel.import_from_pkl(filename='Data/model_8_channels/lower_plenum.pkl')
ex_channel = Channel.import_from_pkl(filename='Data/model_8_channels/external_loop.pkl')
downcomer_channel = Channel.import_from_pkl(filename='Data/model_8_channels/downcomer.pkl')


# channel array
ch_arr_list = [ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8]
ratios = [ch1.mdot_bc,ch2.mdot_bc,ch3.mdot_bc,ch4.mdot_bc,ch5.mdot_bc,ch6.mdot_bc,ch7.mdot_bc,ch8.mdot_bc]

ratios = ratios / sum(ratios) # normalize to 1.0
print("bypass after normalizing first time is ", ratios[-1])
scale_me = 1.0 - PERCENT/100.0
original_sum = sum(ratios[:-1])
ratios[:-1] *= scale_me/original_sum
ratios[-1] = PERCENT/100.0
print("ratios new are: ", ratios)
print("bypass ratio is: ", ratios[-1])
print("sum of ratios is: ", sum(ratios))

ch = ChannelArray(channels=ch_arr_list,
                  coupling_method='ratio_method',
                  flow_ratios=ratios,
                  fluid=fluid,
                  mdot_relaxation=1.0,
                  epsilon=1e-6)
# set bcs for the channel array.
ch.set_bcs(pressure_bc=pressure_bc,
           T_bc=T_bc,
           mdot_bc=mdot_bc,
           tracer_name_value_pairs={}, tracer_bool=False, th_bool=True)

### PRINT CHANNEL RESIDENCE TIMES
print("CHANNEL RESIDENCE TIMES ARE ...")
print(ch.get_channel_residence_time())
print(upper_plenum_channel.get_channel_residence_time())
print(ex_channel.get_channel_residence_time())
print(downcomer_channel.get_channel_residence_time())
print(lower_plenum_channel.get_channel_residence_time())
print("total residence time =", ch.get_channel_residence_time()
      + upper_plenum_channel.get_channel_residence_time()
      + ex_channel.get_channel_residence_time()
      + downcomer_channel.get_channel_residence_time()
      + lower_plenum_channel.get_channel_residence_time())

### INTERFACE SETUPS ###
ch_to_up = ChannelInterface(ch1=ch, ch2=upper_plenum_channel)
up_to_ex = ChannelInterface(ch1=upper_plenum_channel, ch2=ex_channel)
ex_to_downcomer = ChannelInterface(ch1=ex_channel, ch2=downcomer_channel)
downcomer_to_lp = ChannelInterface(ch1=downcomer_channel, ch2=lower_plenum_channel)
lp_to_ch = ChannelInterface(ch1=lower_plenum_channel, ch2=ch)

"""
Steady State Solutions
"""
# st.st. solution for thermal hydraulics
ch.solve_channel_TH(_dt=1e321)
ch_to_up.update_interface_conditions(tracer_bool=False, th_bool=True)

upper_plenum_channel.solve_channel_TH(_dt=1e321)
up_to_ex.update_interface_conditions(tracer_bool=False, th_bool=True)

ex_channel.solve_channel_TH(_dt=1e321)
ex_to_downcomer.update_interface_conditions(tracer_bool=False, th_bool=True)

downcomer_channel.solve_channel_TH(_dt=1e321)
downcomer_to_lp.update_interface_conditions(tracer_bool=False, th_bool=True)

lower_plenum_channel.solve_channel_TH(_dt=1e321)
lp_to_ch.update_interface_conditions(tracer_bool=False, th_bool=True)

# St.St. Solution for tracers
ch_outlet_value = 1.0
iter_num = 0
tracer_to_converge_on = 'c1'
while True:
  iter_num += 1
  ch.solve_all_tracers(_dt=1e321)
  ch_to_up.update_interface_conditions(tracer_bool=True, th_bool=False)

  upper_plenum_channel.solve_all_tracers(_dt=1e321)
  up_to_ex.update_interface_conditions(tracer_bool=True, th_bool=False)

  ex_channel.solve_all_tracers(_dt=1e321)
  ex_to_downcomer.update_interface_conditions(tracer_bool=True, th_bool=False)

  downcomer_channel.solve_all_tracers(_dt=1e321)
  downcomer_to_lp.update_interface_conditions(tracer_bool=True, th_bool=False)

  lower_plenum_channel.solve_all_tracers(_dt=1e321)
  lp_to_ch.update_interface_conditions(tracer_bool=True, th_bool=False)

  _,_,_,_,tracer_outlet_weighted_value = ch.get_outlet_conditions()
  diff = np.abs(1 - tracer_outlet_weighted_value[tracer_to_converge_on]/ch_outlet_value)

  ch_outlet_value = tracer_outlet_weighted_value[tracer_to_converge_on]
  print("Iteration nuimber", iter_num, " and diff is", diff)
  if diff < 1e-12:
    break


### PRINT SOURCES I GUESS ? ###
total_source_times_V = ch.integrate_tracer_source(tracer_name='c1') + upper_plenum_channel.integrate_tracer_source(tracer_name='c1') + lower_plenum_channel.integrate_tracer_source(tracer_name='c1') + downcomer_channel.integrate_tracer_source(tracer_name='c1')

print(ch.integrate_tracer_source(tracer_name='c1')/total_source_times_V)
print(upper_plenum_channel.integrate_tracer_source(tracer_name='c1')/total_source_times_V)
print(lower_plenum_channel.integrate_tracer_source(tracer_name='c1')/total_source_times_V)
print(downcomer_channel.integrate_tracer_source(tracer_name='c1')/total_source_times_V)

"""

  ### TRANSIENTS ###

"""
# TIME SETTINGS
Tstart = -0.1
Tend = 50
nsteps = 2500

# TIMESTEPS
timesteps = np.linspace(Tstart, Tend, nsteps)

# Update old values in channel before starting simulation
ch.update_old_to_most_recent()
upper_plenum_channel.update_old_to_most_recent()
ex_channel.update_old_to_most_recent()
lower_plenum_channel.update_old_to_most_recent()
downcomer_channel.update_old_to_most_recent()

# Transient solver
t_prev = -1e321

# setup beff dict
beta_eff_dict_mc = {}

# Tracer names for stuff
tracer_names = ['c1', 'c2', 'c3', 'c4', 'c5', 'c6']

for t_new in timesteps:
  # DT
  this_dt = t_new - t_prev
  # PRINT TIME INFORMATION
  print("NOW SOLVING AT TIME =", t_new, "| this_dt =", this_dt)
  # CHANGE MASS FLOW RATE AT CHANNEL INLET
  THIS_X_VALUE = t_new
  mdot_bc = MDOT_MAX / 100 * csv_interpolator(csv_df=mdot_csv, x_value=THIS_X_VALUE, x_label='time', y_label='mdot')
  ch.mdot_bc = mdot_bc
  print("\tmdot_bc =", mdot_bc)

  # THERMAL HYDRAULIC INFORMATION
  ch.solve_channel_TH(_dt=this_dt)
  ch_to_up.update_interface_conditions(tracer_bool=False, th_bool=True)

  upper_plenum_channel.solve_channel_TH(_dt=this_dt)
  up_to_ex.update_interface_conditions(tracer_bool=False, th_bool=True)

  ex_channel.solve_channel_TH(_dt=this_dt)
  ex_to_downcomer.update_interface_conditions(tracer_bool=False, th_bool=True)

  downcomer_channel.solve_channel_TH(_dt=this_dt)
  downcomer_to_lp.update_interface_conditions(tracer_bool=False, th_bool=True)

  lower_plenum_channel.solve_channel_TH(_dt=this_dt)
  lp_to_ch.update_interface_conditions(tracer_bool=False, th_bool=True)

  # TRACER TRANSIENT SOLUTION
  ch.solve_all_tracers(_dt=this_dt)
  ch_to_up.update_interface_conditions(tracer_bool=True, th_bool=False)

  upper_plenum_channel.solve_all_tracers(_dt=this_dt)
  up_to_ex.update_interface_conditions(tracer_bool=True, th_bool=False)

  ex_channel.solve_all_tracers(_dt=this_dt)
  ex_to_downcomer.update_interface_conditions(tracer_bool=True, th_bool=False)

  downcomer_channel.solve_all_tracers(_dt=this_dt)
  downcomer_to_lp.update_interface_conditions(tracer_bool=True, th_bool=False)

  lower_plenum_channel.solve_all_tracers(_dt=this_dt)
  lp_to_ch.update_interface_conditions(tracer_bool=True, th_bool=False)

  channel_array_weights = ch.get_channel_tracer_sources(tracer_name='c1')
  up_F = upper_plenum_channel.get_channel_tracer_sources(tracer_name='c1')
  lp_F = lower_plenum_channel.get_channel_tracer_sources(tracer_name='c1')
  downcomer_F = downcomer_channel.get_channel_tracer_sources(tracer_name='c1')

  this_beff_mc, _ = compute_beff_multichannel(channels=ch.channels+[upper_plenum_channel, lower_plenum_channel, downcomer_channel],
                                              weights=channel_array_weights+[up_F, lp_F, downcomer_F],
                                              names=tracer_names)
  beta_eff_dict_mc[t_new] = this_beff_mc

  # Save channel data:
  ch.save_data(_t=t_new)
  upper_plenum_channel.save_data(_t=t_new)
  ex_channel.save_data(_t=t_new)
  lower_plenum_channel.save_data(_t=t_new)
  downcomer_channel.save_data(_t=t_new)

  ch.update_old_to_most_recent()
  upper_plenum_channel.update_old_to_most_recent()
  ex_channel.update_old_to_most_recent()
  lower_plenum_channel.update_old_to_most_recent()
  downcomer_channel.update_old_to_most_recent()

  t_prev = t_new


### MSRE DATA ###
### EXP VALUES ###
msre_integral_worth = 'Data/msre_integral_rod_worth.csv'
msre_data = 'Data/msre_startup_datapoints.csv'

df = pd.read_csv(msre_integral_worth)
data = pd.read_csv(msre_data)

z = df['z'].values
rho = df[' rho'].values
pos = data['pos'].values
ornl_time_data = data['time'].values

values = np.interp(pos, z, rho)
base_insertion = values[0]
startup_base_insertion = values[0]

ornl_delta_rho_data = 1000*(values - base_insertion) # this is the data we calculated based on experimental values

# MAKE PLOT
plt.figure(figsize=(5,3))
plt.plot(ornl_time_data, ornl_delta_rho_data, 'ko', markerfacecolor='w', markersize=3)

### MY VALUES ###
x = np.array(list(beta_eff_dict_mc.keys()))
y_MC = beta_eff_dict_mc[-0.1] - np.array(list(beta_eff_dict_mc.values()))
plt.plot(x, y_MC*10**5, 'k--')
# plt.plot(x[0::50], y[0::50]*10**5, 'ks', markerfacecolor='w')
# plt.plot(x[0::50], y[0::50]*10**5, 'k+', markerfacecolor='w')

# OTHER PLOTTING STUFF
plt.grid()
plt.xlim([-2,50])
plt.ylim([0,300])
print(beta_eff_dict_mc[-0.1])

import pickle as pkl

with open(OUTPUT_FILENAME, "wb") as file:
    pkl.dump(beta_eff_dict_mc, file)
