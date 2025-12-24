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

from Kernels.Helpers import *
from Kernels.NusseltModels import *

from Paraview.Mesher import ParaviewMeshFromChannelArray as PMesh

from Kernels.LumpedCapacitor import *
from Kernels.Helpers import *

"""
Load msre information
"""

# FLUID
fluid = msre_data_dict['fluid']

### BC DATA ###
# mdot_csv = load_csv('Data/msre_flow_coastdown.csv')
mdot_bc = msre_data_dict['mdot_max']
T_bc = msre_data_dict['temp_bc']
pressure_bc = msre_data_dict['pressure_bc']

"""
Load the channels
"""

with open("Data/1152_channels.pkl", "rb") as f:
  obj = pkl.load(f)

# Channel objects
arr: ChannelArray = obj['chArr']
downcomer = obj['DOWN']
external = obj['EX']
upper_p = obj['UP']
lower_p = obj['LP']

# Interfaces
ch_to_up = obj['ch_to_up']
downcomer_to_lp = obj['downcomer_to_lp']
ex_to_downcomer = obj['ex_to_downcomer']
lp_to_ch = obj['lp_to_ch']
up_to_ex = obj['up_to_ex']

"""
Paraview mesh setup
"""
pview = PMesh(arr=arr, filename='thermal_1152.vtu', dx=0.0254 * 2.0/2**0.5, dy=0.0254 * 2.0/2**0.5, writeConductors=True)

"""
Conductor Setup
"""
graphite_k = 80.0
graphite_cp = 1983.0
graphite_multiplier = 1e-5  # W/m3 in fuel -> W/m3 in graphite
# singh et al nonlinear dynamic model results.
salt_core_mass = 1374.0
graphite_core_mass = 3634.0
mass_ratio = graphite_core_mass / salt_core_mass
graphite_node_mass = arr.channels[0].vol_vec[0] * arr.channels[0].rho.T[0] * mass_ratio  # kg of graphite in each node

### Get heat transfer areas
full_ch_perimeter = msre_data_dict['full_ch_perimeter']
half_ch_perimeter = msre_data_dict['half_ch_perimeter']

# LMP
k_LMP = LumpedMaterialProperty(name='thermal conductivity', function_coeffs=[graphite_k])
cp_LMP = LumpedMaterialProperty(name='specific heat', function_coeffs=[graphite_cp])

# Nusselt model
the_nu_model = BasicNusseltModel()

# Attach conductors to all subchannels (not the bypass region though)
for ch in arr.channels:

  if ch.area > 0.036:
    # This is the bypass region
    continue
  if ch.area > 0.00028:
    # This is the full channels
    peri = full_ch_perimeter
    ht_area = ch.L1 / ch.nZones*peri
    gr = LumpedCapacitor(mass=graphite_node_mass,
                        power=-1,
                        h=-1,
                        A=ht_area,
                        C=cp_LMP,
                        thermal_cond=k_LMP,
                        initial_T=900,
                        T_bulk=-1,
                        epsilon=1e-5)

    ch.set_nu_model(nu_model = the_nu_model)
    ChannelConductionBuilder(ch, gr)
  elif ch.area > 0.00014:
    # Half channels
    peri = half_ch_perimeter
    ht_area = ch.L1 / ch.nZones*peri
    gr = LumpedCapacitor(mass=graphite_node_mass,
                         power=-1,
                         h=-1,
                         A=ht_area,
                         C=cp_LMP,
                         thermal_cond=k_LMP,
                         initial_T=900,
                         T_bulk=-1,
                         epsilon=1e-5)
    ChannelConductionBuilder(ch, gr)
    ch.set_nu_model(nu_model = the_nu_model)

  else:
    raise Exception("Channel area unknown.")

# Linking conductors to channels via linkers:
LinkConductorPowersToChannelArray(channelarray=arr, multiplier=1.0,ch_idx_DNI=[1152])
LinkConductorHTCToChannelArray(channelarray=arr, multiplier=1.0,ch_idx_DNI=[1152])
LinkConductorTbulkToChannelArray(channelarray=arr, multiplier=1.0,ch_idx_DNI=[1152])

# Now set the total power:
ChannelArraySetTotalPower(arr=arr, cond_power=8e6*0.06, channel_power=8e6*0.94)

"""
BC setup
"""
# Now set the boundary conditions:
MDOT_BC = msre_data_dict['mdot_max']
PRESSURE_BC = msre_data_dict['pressure_bc']
T_BC = msre_data_dict['temp_bc']

arr.set_bcs(pressure_bc=PRESSURE_BC, T_bc=T_BC, mdot_bc=MDOT_BC, tracer_bool=False, th_bool=True, tracer_name_value_pairs={})
arr.reset_solve()
arr.set_residual_information(energy=1e7, energy_abs=1e-6, print_energy=False)

print(f"MDOT_BC is {MDOT_BC}\nTEMP_BC is f{T_BC}\nPRESSURE_BC is {pressure_bc}")

"""
Solving....
"""
diff = 1e321
n = 0
import time

while n < 10:
  begin_time =  time.time()
  # OLD = copy.deepcopy(arr.temp.T)
  arr.solve_channel_TH(_dt = 1e321, N=31)
  n += 1
  print(f"n = {n}")
  end_time = time.time()
  print(f"time is = {end_time - begin_time}")


