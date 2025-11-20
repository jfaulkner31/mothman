"""
Imports
"""
import numpy as np


"""
Put user inputs here...

Automatically adjusts length of channel to match desired flow tau value (to match residence time in the main part of the core.)

"""
### USER INPUTS ###
TEMP = 900.0
PRESSURE_BC = 101325.0
MSRE_MDOT = 162.0
DOWNCOMER_TAU = 1.3
LP_TAU = 3.8
UP_TAU = 3.9
TOTAL_CIRC_TIME = 25.19
MAIN_CORE_TAU = 9.1725


# FUEL PROPERTIES
CP_FUEL = 1983.0
RHO_0_FUEL = 2715.13
DRHO_DT_FUEL = -0.513
MU_FUEL = 0.00744
K_FUEL = 1.44

# OTHER
GRAVITY = 9.81
TEMP_TOLERANCE = 1.0
MAX_TEMP_ITERATIONS = 10000
NZONES = 50

# DNP DATA (openmc i believe)
TRACER_NAMES = ['c1', 'c2', 'c3', 'c4', 'c5', 'c6']
TRACER_SCHEME = 'upwind'
INITIAL_VALUE_TRACERS = 0.0
DECAY_CONSTS = [0.01334, 0.03273, 0.1208, 0.3029, 0.8498, 2.854]
BETA = [0.000228, 0.001177, 0.001124, 0.002522, 0.001036, 0.000434]
PHI_TRACER = 1.0
RHO_TRACER = 1.0

# GEOMETRY CHANNEL DATA
FULL_CH_AREA = 0.000287587162368000000
FULL_CH_DH = 0.0158276279311347000000
HALF_CH_AREA = 0.000143793581184000000
HALF_CH_DH = 0.0079138139655673500000
FULL_CH_PERIMIETER = 4 * FULL_CH_AREA / FULL_CH_DH
HALF_CH_PERIMETER = 4 * HALF_CH_AREA / HALF_CH_DH


CORE_FULL_FLOW_AREA = 0.36334036056183217 # flow area of the full core.

# OTHER CHANNEL DATA
BYPASS_DH = (71.097 - 70.285)*2.0 / 100

# NUMBER OF CHANNELS IN THE WHOLE CORE
CHANNEL_AREA_RATIO_HORIZONTAL =     np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,0.0,0.0,0.0,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.0,0.0,0.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0],
                                              [0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0],
                                              [0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0],
                                              [0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0],
                                              [0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0],
                                              [0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0],
                                              [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1],
                                              [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1],
                                              [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1],
                                              [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1],
                                              [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1],
                                              [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1],
                                              [0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0],
                                              [0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0],
                                              [0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0],
                                              [0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0],
                                              [0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0],
                                              [0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,0.0,0.0,0.0,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.0,0.0,0.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0],
                                              [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0]])
CHANNEL_AREA_RATIO_VERTICAL = np.transpose(CHANNEL_AREA_RATIO_HORIZONTAL)
NUM_WHOLE_CHANNELS = int(np.sum(CHANNEL_AREA_RATIO_HORIZONTAL == 1.0) * 2.0)
NUM_HALF_CHANNELS = int(np.sum(CHANNEL_AREA_RATIO_HORIZONTAL == 0.5) * 2.0)


"""
Now doing calculations....
"""
print("NOW RUNNING MSRE DATA USER INPUT")
### CALCULATIONS ###
from Subchannel.FluidRelation import FluidRelation
FUEL_FLUID = FluidRelation(cp=CP_FUEL, mu=MU_FUEL, k=K_FUEL, rho_0 = RHO_0_FUEL, drho_dT=DRHO_DT_FUEL)
density = FUEL_FLUID.props_from_P_T(P=101325.0, T=TEMP, prop='rho')

default_length = 1.0

m3s = MSRE_MDOT / density

ch_length = m3s * MAIN_CORE_TAU / CORE_FULL_FLOW_AREA # length of channel
L0 = 0.0
L1 = ch_length
EXTERNAL_TAU = TOTAL_CIRC_TIME - DOWNCOMER_TAU - LP_TAU - UP_TAU - MAIN_CORE_TAU
print("EXTENRAL TAU IS", EXTERNAL_TAU)
print("Channel length L1 is", L1)
print("\n")
# thermal hydraulic flow areas.
LP_AREA = LP_TAU/default_length*m3s
UP_AREA = UP_TAU/default_length*m3s
DOWNCOMER_AREA = DOWNCOMER_TAU/default_length*m3s
EXTERNAL_LOOP_AREA = EXTERNAL_TAU/default_length*m3s
print("LP AREA IS", LP_AREA)
print("UP AREA IS", UP_AREA)
print("DOWNCOMER AREA IS", DOWNCOMER_AREA)
print("EXTERNAL LOOP AREA IS", EXTERNAL_LOOP_AREA)

msre_data_dict = {}

msre_data_dict['tracer_names'] = TRACER_NAMES
msre_data_dict['tracer_scheme'] = TRACER_SCHEME
msre_data_dict['decay_consts'] = DECAY_CONSTS
msre_data_dict['beta'] = BETA
msre_data_dict['phi_tracer'] = PHI_TRACER
msre_data_dict['rho_tracer'] = RHO_TRACER
msre_data_dict['initial_value_tracers'] = INITIAL_VALUE_TRACERS

msre_data_dict['full_ch_area'] = FULL_CH_AREA
msre_data_dict['full_ch_dh'] = FULL_CH_DH
msre_data_dict['half_ch_area'] = HALF_CH_AREA
msre_data_dict['half_ch_dh'] = HALF_CH_DH
msre_data_dict['full_ch_perimeter'] = FULL_CH_PERIMIETER
msre_data_dict['half_ch_perimeter'] = HALF_CH_PERIMETER


msre_data_dict['core_full_flow_area'] = CORE_FULL_FLOW_AREA

msre_data_dict['temp_bc'] = TEMP
msre_data_dict['pressure_bc'] = PRESSURE_BC
msre_data_dict['mdot_max'] = MSRE_MDOT

msre_data_dict['lp_area'] = LP_AREA
msre_data_dict['up_area'] = UP_AREA
msre_data_dict['dc_area'] = DOWNCOMER_AREA
msre_data_dict['ex_area'] = EXTERNAL_LOOP_AREA
msre_data_dict['bypass_dh'] = BYPASS_DH
msre_data_dict['L0'] = L0
msre_data_dict['L1'] = L1

msre_data_dict['fluid'] = FUEL_FLUID

msre_data_dict['gravity'] = GRAVITY
msre_data_dict['temp_tol'] = TEMP_TOLERANCE
msre_data_dict['max_temp_iterations'] = MAX_TEMP_ITERATIONS
msre_data_dict['nZones'] = NZONES

msre_data_dict['channel_area_ratio_horizontal'] = CHANNEL_AREA_RATIO_HORIZONTAL
msre_data_dict['channel_area_ratio_vertical'] = CHANNEL_AREA_RATIO_VERTICAL

msre_data_dict['num_whole_channels'] = NUM_WHOLE_CHANNELS
msre_data_dict['num_half_channels'] = NUM_HALF_CHANNELS


