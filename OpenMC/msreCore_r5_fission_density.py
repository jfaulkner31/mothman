############################################
### TODO BLOCK
############################################
# tally over 252g difusion coeffs and then condense.
# also tally control rod regions and fuel regions seperately..?
# correct lower part of control rod to extend all the way to the bottom of the graphite - otherwise introduces needless heterogeneity.
# change meshing - have mesh start where control rod starts..... oops
# 252g transport xs tallying.
# tally fuel salt outside of each control rod?

# same radial tallying loadout but maybe less axial zones.
# 1 zone in crod layer where resting rods are
# another zone in crod layer where inserted rod sits.
# then 6 equidistant zones until bottom of core.

# rev4 update:
# added a few more channels on periphary
# changed channel dimensions to thermally expanded dimensions in MSRE benchmark
# changed core can, graphite lattice, etc. dimensions to match th. expanded dimensions in MSRE benchmark
#

############################################
### HEADER
############################################
print('####################################################################################################################################')
print('####################################################################################################################################')
print('####################################################################################################################################')
print('####################################################################################################################################')
############################################
### IMPORTS
############################################
# %matplotlib inline
import openmc
import numpy as np
import matplotlib.pyplot as plt
import openmc.mgxs as mgxs
from openmcToGriffinMakeTallies import *
import sys

# openmc.config['cross_sections'] = "/mnt/c/openMC_code/endfb-vii.1-hdf5/cross_sections.xml"
# openmc.config['cross_sections'] = "/mnt/n/endfb-vii.1-hdf5/cross_sections.xml"

############################################
### SYS ARGUMENTS
############################################

checkArgsLen = len(sys.argv)

if checkArgsLen != 15:
  print("sysArgs[:] = ", sys.argv[:])
  print("length of args = ", checkArgsLen)
  print('rod 1 depth inches')
  print('rod 2 depth inchES')
  print('NPG')
  print('NSK')
  print('NACTIVE')
  print('NTHREADS')
  print('233 or 235 based on fuel type')
  print('fuel temp')
  print('other temps')
  print('graphite temp')
  print("enter '1' for XS generation ON and '0' for XS generation OFF")
  print("enter group structure index from list in xs gen section")
  # print("number of axial intervals for xs generation over defined mesh")
  raise Exception("Wrong size of input argv's")
rod_1_user_input_depth_inches = float(sys.argv[1]) # (cm measured from top of graphite) aka insertion depth ---> ctrl+f for precise usage of these.
rod_2_user_input_depth_inches = float(sys.argv[2]) # (cm measured from top of graphite) aka insertion depth
npg =  int(sys.argv[3])
nsk =  int(sys.argv[4])
nact = int(sys.argv[5])
nthreads = int(sys.argv[6])
fuel_type = int(sys.argv[7])
fuel_temperature = float(sys.argv[8]) # temperature flibe
temperature = float(sys.argv[9]) # temperature all other materials
graphite_temperature =float(sys.argv[10]) # temperature graphite
fuel_density_type = float(sys.argv[11]) # type 1 for constant and type 2 for variable
xs_gen_logical = int(sys.argv[12])
user_input_group_struct = int(sys.argv[13]) # group structure choice from list
xs_gen_axial_intervals = int(sys.argv[14])

print()
print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
print("User input to insert rod 1 to ", rod_1_user_input_depth_inches, "inches.")
print("total rod 1 depth (below graphite top) is 2.75 + ", rod_1_user_input_depth_inches, " inches")
print("User input to insert rod 2 to ", rod_2_user_input_depth_inches, "inches.")
print("total rod 2 depth (below graphite top) is 2.75 + ", rod_2_user_input_depth_inches, " inches")
print("User Input NPG: ", npg)
print("User Input NSK: ", nsk)
print("User Input NACT: ", nact)
print("Number of threads: ", nthreads)
print("Fuel Type - U"+str(fuel_type))
print("Fuel temperature - ", fuel_temperature)
print("Other mat temperatures (besides thermal shield and insulation) - ", temperature)
print("Graphite Temperature - ", graphite_temperature)
print("fuel density type: 1 if constant 2.3275 OR 2 if based on formulation.", fuel_density_type)
print("XS Generation (1 for on and 0 for off ) = ", xs_gen_logical)
print("Energy group structure choice = ", user_input_group_struct)
print("Number xs axial intervals", xs_gen_axial_intervals)
print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
print()
print()
############################################
### USER PARAMS
############################################
thermal_shield_temperature = 305  # from benchmark pdf page 75 3.1.4 NEA
insulation_temperature = 305 # from benchmark pdf page 75 3.1.4 NEA
if fuel_density_type == 2:
  # fuel density from Zanetti et al. PHYSOR paper - https://thoriumenergyalliance.com/wp-content/uploads/2020/02/Dynamics-Simulations-of-MSRE.pdf
  fuel_density = 2.575-5.13e-4*(float(fuel_temperature)-273.15)# g/cm3  rho0 - a*T(C) - if T is in kelvin, we subtract 273.15 -
  other_fuel_density = 2.575-5.13e-4*(float(temperature)-273.15) # fuel that is out of core

elif fuel_density_type == 1:
  fuel_density = 2.3275 # msre benchmark desc
  other_fuel_density = 2.3275 # fuel that is in-core

else:
  raise Exception("unknown fuel density type")

graphite_density = 1.86 # from msre benchmark
helium_density = 0.000072 * (911/temperature) # density at 911K from msre benchmark
poison_density = 5.873 # benchmark room temp density
inor_density = 8.7745 # from msre benchmark
ss316_density = 8.03 # from scale ss316 and from msre benchmark
inconel_density = 8.192 #from msre benchmark
tshield_density = 4.42 # from msre benchmark
insulation_density = 0.112 # msre benchmark again
cellgas_density = 0.000505 * (911/temperature) # from benchmark description density at 911 K
cellgas_density_in_core = 0.000505 * (911/graphite_temperature)
print('Fuel Density Is: ', fuel_density)

############################################
### MATERIALS
############################################
#### MSRE MATS FROM SPREADSHEET
#### SEE rev2 of origen spreadsheet for calculation -> Tech -> senior design tutorials -> origen -> rev2 spreadsheet
#### made based on sample FP-17-1 from ORNL-2076 report (Chemical Aspects of MSRE Operations. Roy E. Thoma)
#### note that this had to be slightly corrected since sum(x) was greater than 100 wt% so we reduced the F-19 isotope to make it match up
#### representative of U-233 operations, need to use for getting dnp data for griffin simulations
#### use same density relation as otherwise.

# TODO MSRE FROM OPENMSRE
# Mat comps stuff calculated based on values from rhanema and students thesis - https://repository.gatech.edu/server/api/core/bitstreams/06739ab4-be94-4ec0-85c1-03aa1f56d333/content
# f1=openmc.Material(name='f1')
# f1.add_nuclide('Li6',1.31480070E-05)
# f1.add_nuclide('Li7', 0.262960140146177)
# f1.add_nuclide('Be9',1.1863E-01)
# f1.add_nuclide('Zr90',1.0543E-02)
# f1.add_nuclide('Zr91',2.2991E-03)
# f1.add_nuclide('Zr92',3.5142E-03)
# f1.add_nuclide('Zr94',3.5613E-03)
# f1.add_nuclide('Zr96',5.7375E-04)
# f1.add_nuclide('Hf174',8.3786E-10)
# f1.add_nuclide('Hf176',2.7545E-08)
# f1.add_nuclide('Hf177',9.7401E-08)
# f1.add_nuclide('Hf178',1.4285E-07)
# f1.add_nuclide('Hf179',7.1323E-08)
# f1.add_nuclide('Hf180',1.8370E-07)
# f1.add_nuclide('U234',1.034276246E-05)
# f1.add_nuclide('U235',1.009695816E-03)
# f1.add_nuclide('U236',4.227809892E-06)
# f1.add_nuclide('U238',2.168267822E-03)
# f1.add_nuclide('Fe54',2.8551E-06)
# f1.add_nuclide('Fe56',4.4818E-05)
# f1.add_nuclide('Fe57',1.0350E-06)
# f1.add_nuclide('Fe58',1.3775E-07)
# f1.add_nuclide('Cr50',2.1224E-06)
# f1.add_nuclide('Cr52',4.0928E-05)
# f1.add_nuclide('Cr53',4.6409E-06)
# f1.add_nuclide('Cr54',1.1552E-06)
# f1.add_nuclide('Ni58',5.8597E-06)
# f1.add_nuclide('Ni60',2.2571E-06)
# f1.add_nuclide('Ni61',9.8117E-08)
# f1.add_nuclide('Ni62',3.1284E-07)
# f1.add_nuclide('Ni64',7.9671E-08)
# f1.add_nuclide('O16',5.1437E-04)
# f1.add_nuclide('O17',1.8927E-07+9.6440E-07)
# f1.add_nuclide('F19',5.9409E-01)
# f1.set_density('g/cm3', fuel_density)
# f1.temperature = temperature
# f1.id = 1
### MSRE TODO FROM BENCHMARK DESCRIPTION u235 type
if fuel_type == 235:
  f1=openmc.Material(name='f1') # f1 is in-core salt
  f1.add_nuclide('Li6', 1.0944E-06, 'ao')
  f1.add_nuclide('Li7', 2.1888E-02, 'ao')
  f1.add_element('Be', 9.8747E-03, 'ao')
  f1.add_element('Zr', 1.7056E-03, 'ao')
  f1.add_element('Hf', 4.3588E-08, 'ao')
  f1.add_nuclide('U234', 8.6090E-07, 'ao')
  f1.add_nuclide('U235', 8.4044E-05, 'ao')
  f1.add_nuclide('U236', 3.5191E-07, 'ao')
  f1.add_nuclide('U238', 1.8048E-04, 'ao')
  f1.add_element('Fe', 4.0658E-06, 'ao')
  f1.add_element('Cr', 7.5479E-07, 'ao')
  f1.add_element('Ni', 7.1646E-07, 'ao')
  f1.add_element('O', 4.2927E-05, 'ao')
  f1.add_element('F', 4.9450E-02, 'ao')
  f1.set_density('g/cm3', fuel_density)
  f1.temperature = fuel_temperature
  f1.id = 1
  # f2 is fuel in upper/lower plena, downcomer, upcomer
  f2 = openmc.Material(name='f2')
  f2.add_nuclide('Li6', 1.0944E-06, 'ao')
  f2.add_nuclide('Li7', 2.1888E-02, 'ao')
  f2.add_element('Be', 9.8747E-03, 'ao')
  f2.add_element('Zr', 1.7056E-03, 'ao')
  f2.add_element('Hf', 4.3588E-08, 'ao')
  f2.add_nuclide('U234', 8.6090E-07, 'ao')
  f2.add_nuclide('U235', 8.4044E-05, 'ao')
  f2.add_nuclide('U236', 3.5191E-07, 'ao')
  f2.add_nuclide('U238', 1.8048E-04, 'ao')
  f2.add_element('Fe', 4.0658E-06, 'ao')
  f2.add_element('Cr', 7.5479E-07, 'ao')
  f2.add_element('Ni', 7.1646E-07, 'ao')
  f2.add_element('O', 4.2927E-05, 'ao')
  f2.add_element('F', 4.9450E-02, 'ao')
  f2.set_density('g/cm3', other_fuel_density)
  f2.temperature = temperature
  f2.id = 2
elif fuel_type == 233:
# MSRE TODO U233 TYPE MSR FROM BURKE AND RHANEMA BASED ON SPREADSHEET GENERATION VALUES...
  f1=openmc.Material(name='f1')
  f1.add_nuclide('Li6', 1.52823674612924E-05, 'ao')
  f1.add_nuclide('Li7', 0.262030707719481, 'ao')
  f1.add_nuclide('F19', 0.593727147152027, 'ao')
  f1.add_nuclide('Be9', 0.122613146989518, 'ao')
  f1.add_nuclide('Zr90', 0.0108485211668156, 'ao')
  f1.add_nuclide('Zr91', 0.00236579995124726, 'ao')
  f1.add_nuclide('Zr92', 0.00361617372227188, 'ao')
  f1.add_nuclide('Zr94', 0.00366467051271634, 'ao')
  f1.add_nuclide('Zr96', 0.000590395709758674, 'ao')
  f1.add_nuclide('U235', 1.26453979317652E-05, 'ao')
  f1.add_nuclide('U238', 3.18617710794216E-05, 'ao')
  f1.add_nuclide('U234', 3.69921530780425E-05, 'ao')
  f1.add_nuclide('U236', 3.65735041963589E-07, 'ao')
  f1.add_nuclide('U233', 0.000446289651571172, 'ao')
  f1.set_density('g/cm3', fuel_density)
  f1.temperature = fuel_temperature
  f1.id = 1
  # f2 is fuel in upper/lower plena, downcomer, upcomer
  f2=openmc.Material(name='f2')
  f2.add_nuclide('Li6', 1.52823674612924E-05, 'ao')
  f2.add_nuclide('Li7', 0.262030707719481, 'ao')
  f2.add_nuclide('F19', 0.593727147152027, 'ao')
  f2.add_nuclide('Be9', 0.122613146989518, 'ao')
  f2.add_nuclide('Zr90', 0.0108485211668156, 'ao')
  f2.add_nuclide('Zr91', 0.00236579995124726, 'ao')
  f2.add_nuclide('Zr92', 0.00361617372227188, 'ao')
  f2.add_nuclide('Zr94', 0.00366467051271634, 'ao')
  f2.add_nuclide('Zr96', 0.000590395709758674, 'ao')
  f2.add_nuclide('U235', 1.26453979317652E-05, 'ao')
  f2.add_nuclide('U238', 3.18617710794216E-05, 'ao')
  f2.add_nuclide('U234', 3.69921530780425E-05, 'ao')
  f2.add_nuclide('U236', 3.65735041963589E-07, 'ao')
  f2.add_nuclide('U233', 0.000446289651571172, 'ao')
  f2.set_density('g/cm3', other_fuel_density)
  f2.temperature = temperature
  f2.id = 2
else:
  throwErrorString = "Fuel type is wrong type - expected either numerical 233 or 235 entries but instead got value: "+str(fuel_type)
  raise Exception(throwErrorString)




# graphite borated from benchmark description
g1= openmc.Material(name='g1')
g1.add_nuclide('C0', 9.2789E-02, 'ao')
g1.add_nuclide('B10', 1.6412E-08,'ao')
g1.add_nuclide('B11', 6.6060E-08,'ao')
g1.add_element('V', 1.9690E-07,'ao')
g1.add_element('S', 1.7378E-07,'ao')
g1.add_element('O', 1.7235E-07,'ao')
g1.add_element('Si', 5.3518E-08,'ao')
g1.add_element('Al', 3.7225E-08,'ao')
g1.add_element('Fe', 3.6679E-09,'ao')
g1.add_element('Ti', 1.2688E-09,'ao')
g1.add_element('Mg', 1.0027E-09,'ao')
g1.add_element('Ca', 4.3868E-10,'ao')
g1.add_s_alpha_beta('c_Graphite', 1.0)
g1.set_density("g/cm3", graphite_density)
g1.temperature = graphite_temperature
g1.id=999

### POISON from benchmark
poison = openmc.Material(name='poison')
poison.add_element('Gd', 1.3310E-02,'ao')
poison.add_element('Al', 2.0280E-02,'ao')
poison.add_element('O', 5.0384E-02,'ao')
poison.set_density('g/cm3', poison_density)
poison.temperature = temperature
poison.id = 10000

# currently from scale wth 8.3g/cc
# from benchmark description
inconel = openmc.Material(name="inconel_in_core")
inconel.add_nuclide('C0', 3.2017E-04, 'ao')
inconel.add_element('Mn', 3.0625E-04, 'ao')
inconel.add_element('P', 2.3279E-05, 'ao')
inconel.add_element('S', 2.2486E-05, 'ao')
inconel.add_element('Si', 5.9905E-04, 'ao')
inconel.add_element('Cr', 1.7565E-02, 'ao')
inconel.add_element('Ni', 4.2998E-02, 'ao')
# inconel.add_element('Mo', 1.5283E-03, 'ao') # remove Mo and boron for simplicity in homogenizing
inconel.add_element('Nb', 2.6517E-03, 'ao')
inconel.add_element('Ti', 9.0383E-04, 'ao')
inconel.add_element('Al', 8.9080E-04, 'ao')
inconel.add_element('Co', 8.1568E-04, 'ao')
#inconel.add_nuclide('B10', 5.3090E-06, 'ao') # remove Mo and boron for simplicity in homogenizing
#inconel.add_nuclide('B11', 2.1369E-04, 'ao') # remove Mo and boron for simplicity in homogenizing
inconel.add_element('Cu', 2.2694E-04, 'ao')
inconel.add_element('Ta', 1.3283E-05, 'ao')
inconel.add_element('Fe', 1.4426E-02, 'ao')
inconel.temperature=temperature
inconel.set_density('g/cm3', inconel_density)
inconel.id=10001

# inconel that is in the core - such as in the control rods. same as graphite temperature
inconel_in_core = openmc.Material(name="inconel_in_core")
inconel_in_core.add_nuclide('C0', 3.2017E-04, 'ao')
inconel_in_core.add_element('Mn', 3.0625E-04, 'ao')
inconel_in_core.add_element('P', 2.3279E-05, 'ao')
inconel_in_core.add_element('S', 2.2486E-05, 'ao')
inconel_in_core.add_element('Si', 5.9905E-04, 'ao')
inconel_in_core.add_element('Cr', 1.7565E-02, 'ao')
inconel_in_core.add_element('Ni', 4.2998E-02, 'ao')
# inconel.add_element('Mo', 1.5283E-03, 'ao') # remove Mo and boron for simplicity in homogenizing
inconel_in_core.add_element('Nb', 2.6517E-03, 'ao')
inconel_in_core.add_element('Ti', 9.0383E-04, 'ao')
inconel_in_core.add_element('Al', 8.9080E-04, 'ao')
inconel_in_core.add_element('Co', 8.1568E-04, 'ao')
#inconel.add_nuclide('B10', 5.3090E-06, 'ao') # remove Mo and boron for simplicity in homogenizing
#inconel.add_nuclide('B11', 2.1369E-04, 'ao') # remove Mo and boron for simplicity in homogenizing
inconel_in_core.add_element('Cu', 2.2694E-04, 'ao')
inconel_in_core.add_element('Ta', 1.3283E-05, 'ao')
inconel_in_core.add_element('Fe', 1.4426E-02, 'ao')
inconel_in_core.temperature=graphite_temperature
inconel_in_core.set_density('g/cm3', inconel_density)
inconel_in_core.id=20001

# from benchmark description
ss316 = openmc.Material(name="ss316")
ss316.add_nuclide("C0", 3.1384E-04, "ao")
ss316.add_element('Mc', 1.7154E-03, 'ao')
ss316.add_element('P', 6.8457E-05, 'ao')
ss316.add_element('S', 4.4083E-05, 'ao')
ss316.add_element('Si', 1.2583E-03, 'ao')
ss316.add_element('Cr', 1.7218E-02, 'ao')
ss316.add_element('Ni', 8.0281E-03, 'ao')
ss316.add_element('N', 3.3641E-04, 'ao')
ss316.add_element('Fe', 5.7371E-02, 'ao')
ss316.temperature=temperature
ss316.set_density('g/cm3', ss316_density)
ss316.id = 10002

ss316_in_core = openmc.Material(name="ss316_in_core")
ss316_in_core.add_nuclide("C0", 3.1384E-04, "ao")
ss316_in_core.add_element('Mc', 1.7154E-03, 'ao')
ss316_in_core.add_element('P', 6.8457E-05, 'ao')
ss316_in_core.add_element('S', 4.4083E-05, 'ao')
ss316_in_core.add_element('Si', 1.2583E-03, 'ao')
ss316_in_core.add_element('Cr', 1.7218E-02, 'ao')
ss316_in_core.add_element('Ni', 8.0281E-03, 'ao')
ss316_in_core.add_element('N', 3.3641E-04, 'ao')
ss316_in_core.add_element('Fe', 5.7371E-02, 'ao')
ss316_in_core.temperature=graphite_temperature
ss316_in_core.set_density('g/cm3', ss316_density)
ss316_in_core.id = 20002

# cellgas density
helium = openmc.Material(name="helium")
helium.add_nuclide('He4', 1, 'ao')
helium.temperature=temperature
helium.set_density('g/cm3', helium_density)
helium.id = 10003

# cell gas from benchmark
cellgas=openmc.Material(name='cellgas')
cellgas.add_nuclide('N14', 2.0627e-5, 'ao')
cellgas.add_nuclide('O16', 9.5040e-7, 'ao')
cellgas.temperature=temperature
cellgas.set_density('g/cm3', cellgas_density)
cellgas.id= 10015

# cell gas from benchmark
cellgas_in_core=openmc.Material(name='cellgas_in_core')
cellgas_in_core.add_nuclide('N14', 2.0627e-5, 'ao')
cellgas_in_core.add_nuclide('O16', 9.5040e-7, 'ao')
cellgas_in_core.temperature=graphite_temperature
cellgas_in_core.set_density('g/cm3', cellgas_density_in_core)
cellgas_in_core.id= 20015

# inor from benchmark
inor = openmc.Material(name="inor")
inor.add_element('Ni', 5.9653E-02,'ao')
inor.add_element('Mo', 9.1243E-03,'ao')
inor.add_element('Cr', 6.9316E-03,'ao')
inor.add_element('Fe', 4.6099E-03,'ao')
inor.add_nuclide('C0', 2.5721E-04,'ao')
inor.add_element('Ti', 2.1992E-04,'ao')
inor.add_element('Al', 3.9015E-04,'ao')
inor.add_element('S', 2.6263E-05,'ao')
inor.add_element('Mn', 7.6645E-04,'ao')
inor.add_element('Si', 1.4993E-03,'ao')
inor.add_element('Cu', 2.3192E-04,'ao')
inor.add_nuclide('B10', 7.7508E-06,'ao')
inor.add_nuclide('B11', 3.1198E-05,'ao')
inor.add_element('W', 1.1452E-04,'ao')
inor.add_element('P', 2.0392E-05,'ao')
inor.add_element('Co', 1.4290E-04,'ao')
inor.set_density('g/cm3', inor_density)
inor.temperature = temperature
inor.id = 10004

# inor from benchmark
inor_in_core = openmc.Material(name="inor_in_core")
inor_in_core.add_element('Ni', 5.9653E-02,'ao')
inor_in_core.add_element('Mo', 9.1243E-03,'ao')
inor_in_core.add_element('Cr', 6.9316E-03,'ao')
inor_in_core.add_element('Fe', 4.6099E-03,'ao')
inor_in_core.add_nuclide('C0', 2.5721E-04,'ao')
inor_in_core.add_element('Ti', 2.1992E-04,'ao')
inor_in_core.add_element('Al', 3.9015E-04,'ao')
inor_in_core.add_element('S', 2.6263E-05,'ao')
inor_in_core.add_element('Mn', 7.6645E-04,'ao')
inor_in_core.add_element('Si', 1.4993E-03,'ao')
inor_in_core.add_element('Cu', 2.3192E-04,'ao')
inor_in_core.add_nuclide('B10', 7.7508E-06,'ao')
inor_in_core.add_nuclide('B11', 3.1198E-05,'ao')
inor_in_core.add_element('W', 1.1452E-04,'ao')
inor_in_core.add_element('P', 2.0392E-05,'ao')
inor_in_core.add_element('Co', 1.4290E-04,'ao')
inor_in_core.set_density('g/cm3', inor_density)
inor_in_core.temperature = graphite_temperature
inor_in_core.id = 20004

# mixed lower plenum ---
mixedLowerPlenaMat = openmc.Material.mix_materials([inor,f2], [9.20/100, 90.8/100], 'vo')
mixedLowerPlenaMat.temperature = temperature
mixedLowerPlenaMat.id = 10007
mixedLowerPlenaMat.name = 'mixedLowerPlenaMat'

# thermal shield
tshield = openmc.Material(name='tshield')
tshield.add_element('Fe'  , 4.1902E-02 ,'ao')
tshield.add_element('K'   , 1.3995E-05, 'ao')
tshield.add_nuclide('C0'  , 1.9679E-03, 'ao')
tshield.add_nuclide('B10', 3.3854E-07, 'ao')
tshield.add_nuclide('B11', 1.3627E-06, 'ao')
tshield.add_nuclide('N14', 8.8896E-06, 'ao')
tshield.add_nuclide('H1' , 3.3248E-02, 'ao')
tshield.add_nuclide('O16', 1.6662E-02, 'ao')
tshield.temperature = thermal_shield_temperature
tshield.set_density('g/cm3', tshield_density)
tshield.id = 10009

# insulation material
insulation_mat = openmc.Material(name='insulation')
insulation_mat.add_element('Ca', 2.9373E-05, 'ao')
insulation_mat.add_element('Fe', 1.2589E-04, 'ao')
insulation_mat.add_nuclide('Al27' ,1.2589E-04, 'ao')
insulation_mat.add_element('Si', 3.3570E-04, 'ao')
insulation_mat.add_nuclide('O16'  ,1.3428E-03, 'ao')
insulation_mat.add_nuclide('H1'   ,8.3924E-04, 'ao')
insulation_mat.temperature = insulation_temperature
# insulation_mat.set_density('g/cm3', insulation_density) ### TODO turn off density setting - let it be decided based on benchmark desc.
insulation_mat.id = 10010

materials=openmc.Materials((f1, f2, g1, inconel, poison, helium, ss316, inor, mixedLowerPlenaMat, tshield, insulation_mat, cellgas,
                            inor_in_core, ss316_in_core, cellgas_in_core))
materials.export_to_xml()



############################################
### GRAPHITE FUEL CELL
############################################
# some definitions
rad = 0.50885 # 0.508  cm cold dimension
L = 5.08339    # 5.08 cm cold dimension
W = 3.05309    # 3.048 cm cold dimension
rec_len = W - 2*rad
d = (L-W)/2 + rad
d2 = L -d
# 8 cylinders that go around the unit cell
c1 = openmc.ZCylinder(x0=d-L/2,y0=L-L/2,r=rad)
c2 = openmc.ZCylinder(x0=d2-L/2,y0=L-L/2,r=rad)
c3 = openmc.ZCylinder(x0=L-L/2,y0=d2-L/2,r=rad)
c4 = openmc.ZCylinder(x0=L-L/2,y0=d-L/2,r=rad)
c5 = openmc.ZCylinder(x0=d2-L/2,y0=0-L/2,r=rad)
c6 = openmc.ZCylinder(x0=d-L/2,y0=0-L/2,r=rad)
c7 = openmc.ZCylinder(x0=0-L/2,y0=d-L/2,r=rad)
c8 = openmc.ZCylinder(x0=0-L/2,y0=d2-L/2,r=rad)
# planes at the x=0 positions (shifted by L/2)
p1 = openmc.XPlane(x0=0-L/2)
p2 = openmc.YPlane(y0=L-L/2)
p3 = openmc.XPlane(x0=L-L/2)
p4 = openmc.YPlane(y0=0-L/2)
# planes at the y=0 positions shifted by L/2
p5 = openmc.XPlane(x0=rad  -L/2)
p6 = openmc.YPlane(y0=L-rad  -L/2)
p7 = openmc.XPlane(x0=L-rad -L/2)
p8 = openmc.YPlane(y0=rad -L/2)
# planes seperating the pure rectangle region and the circular regions shifted by L/2
p9 = openmc.XPlane(x0=d-L/2)
p10 = openmc.XPlane(x0=d2-L/2)
p11 = openmc.YPlane(y0=d2-L/2)
p12 = openmc.YPlane(y0=d-L/2)


inside_channel_1 = ( (-c1 | -c2) & (+p6 & -p2) ) | (+p9 & -p10 & +p6 & -p2)  # TOP CHANNEL
inside_channel_2 = ( (-c3 | -c4) & (+p7 & -p3) ) | (+p12 & -p11 & +p7 & -p3) # RIGHT CHANNEL
inside_channel_3 = ( (-c5 | -c6) & (+p4 & -p8) ) | (+p9 & -p10 & +p4 & -p8)  # BOTTOM
inside_channel_4 = ( (-c7 | -c8) & (+p1 & -p5) ) | (+p12 & -p11 & +p1 & -p5) # LEFT
graphite_bounds = ~inside_channel_1 & ~inside_channel_2 & ~inside_channel_3 & ~inside_channel_4


# make fuel cell filled with material f1
fuel_cell = openmc.Cell(fill=f1, region=(    (inside_channel_1 | inside_channel_2 | inside_channel_3 | inside_channel_4) ))
graphite_channel = openmc.Cell(fill=g1,region=(~inside_channel_1 & ~inside_channel_2 & ~inside_channel_3 & ~inside_channel_4))

u = openmc.Universe(cells=[fuel_cell,graphite_channel])

graphite_cell = openmc.Cell(fill=g1)
G = openmc.Universe(cells=(graphite_cell,))


# u.plot(width=(13,13), color_by='material', colors={f1:'red',g1:'pink'})
# G.plot(width=(10,10))

############################################
### SPECIAL FUEL-GRAPHITE CELLS
############################################
# TOP RIGHT QUARTER OF CORE SPECIAL CELLS AND UNIVERSES
# these are color coded based on Figure 3.6 of the MSRE benchmark to match the channel geometry! https://www.osti.gov/servlets/purl/1617123

blue_fuel_cell_TR = openmc.Cell(fill=f1, region=(    (inside_channel_2 | inside_channel_3 | inside_channel_4) ))
blue_grap_cell_TR = openmc.Cell(fill=g1,region=(~inside_channel_2 & ~inside_channel_3 & ~inside_channel_4 | inside_channel_1))
blue_TR = openmc.Universe(cells=[blue_fuel_cell_TR, blue_grap_cell_TR], name='blue_TR')

purple_fuel_cell_TR = openmc.Cell(fill=f1, region=(    (inside_channel_3 | inside_channel_4) ))
purple_grap_cell_TR = openmc.Cell(fill=g1,region=(  ~inside_channel_3 & ~inside_channel_4  | inside_channel_1 | inside_channel_2  ))
purple_TR = openmc.Universe(cells=[purple_fuel_cell_TR, purple_grap_cell_TR], name='purple_TR')

green_fuel_cell_TR = openmc.Cell(fill=f1, region=(    (inside_channel_3) ))
green_grap_cell_TR = openmc.Cell(fill=g1,region=(~inside_channel_3 | inside_channel_1 | inside_channel_2 | inside_channel_4))
green_TR = openmc.Universe(cells=[green_fuel_cell_TR, green_grap_cell_TR], name='green_TR')

yellow_fuel_cell_TR = openmc.Cell(fill=f1, region=(    (inside_channel_4) ))
yellow_grap_cell_TR = openmc.Cell(fill=g1,region=(~inside_channel_4 | inside_channel_1 | inside_channel_2 | inside_channel_3))
yellow_TR = openmc.Universe(cells=[yellow_fuel_cell_TR, yellow_grap_cell_TR], name='yellow_TR')

brown_fuel_cell_TR = openmc.Cell(fill=f1, region=(    (inside_channel_1 | inside_channel_3 | inside_channel_4) ))
brown_grap_cell_TR = openmc.Cell(fill=g1,region=(~inside_channel_1 & ~inside_channel_3 & ~inside_channel_4 | inside_channel_2))
brown_TR = openmc.Universe(cells=[brown_fuel_cell_TR, brown_grap_cell_TR], name='brown_TR')

### BOTTOM RIGHT QUARTER OF CORE
blue_fuel_cell_BR = openmc.Cell(fill=f1, region=(    (inside_channel_2 | inside_channel_1 | inside_channel_4) ))
blue_grap_cell_BR = openmc.Cell(fill=g1,region=(~inside_channel_2 & ~inside_channel_1 & ~inside_channel_4 | inside_channel_3))
blue_BR = openmc.Universe(cells=[blue_fuel_cell_BR, blue_grap_cell_BR], name='blue_BR')

purple_fuel_cell_BR = openmc.Cell(fill=f1, region=(    (inside_channel_1 | inside_channel_4) ))
purple_grap_cell_BR = openmc.Cell(fill=g1,region=(  ~inside_channel_1 & ~inside_channel_4  | inside_channel_3 | inside_channel_2  ))
purple_BR = openmc.Universe(cells=[purple_fuel_cell_BR, purple_grap_cell_BR], name='purple_BR')

green_fuel_cell_BR = openmc.Cell(fill=f1, region=(    (inside_channel_1) ))
green_grap_cell_BR = openmc.Cell(fill=g1,region=(~inside_channel_1 | inside_channel_3 | inside_channel_2 | inside_channel_4))
green_BR = openmc.Universe(cells=[green_fuel_cell_BR, green_grap_cell_BR], name='green_BR')

yellow_fuel_cell_BR = openmc.Cell(fill=f1, region=(    (inside_channel_4) ))
yellow_grap_cell_BR = openmc.Cell(fill=g1,region=(~inside_channel_4 | inside_channel_1 | inside_channel_2 | inside_channel_3))
yellow_BR = openmc.Universe(cells=[yellow_fuel_cell_BR, yellow_grap_cell_BR], name='yellow_BR')

brown_fuel_cell_BR = openmc.Cell(fill=f1, region=(    (inside_channel_1 | inside_channel_3 | inside_channel_4) ))
brown_grap_cell_BR = openmc.Cell(fill=g1,region=(~inside_channel_1 & ~inside_channel_3 & ~inside_channel_4 | inside_channel_2))
brown_BR = openmc.Universe(cells=[brown_fuel_cell_BR, brown_grap_cell_BR], name='brown_BR')


### BOTTOM LEFT QUARTER OF CORE
blue_fuel_cell_BL = openmc.Cell(fill=f1, region=(    (inside_channel_2 | inside_channel_1 | inside_channel_4) ))
blue_grap_cell_BL = openmc.Cell(fill=g1,region=(~inside_channel_2 & ~inside_channel_1 & ~inside_channel_4 | inside_channel_3))
blue_BL = openmc.Universe(cells=[blue_fuel_cell_BL, blue_grap_cell_BL], name='blue_BL')

purple_fuel_cell_BL = openmc.Cell(fill=f1, region=(    (inside_channel_1 | inside_channel_2) ))
purple_grap_cell_BL = openmc.Cell(fill=g1,region=(  ~inside_channel_1 & ~inside_channel_2  | inside_channel_3 | inside_channel_4  ))
purple_BL = openmc.Universe(cells=[purple_fuel_cell_BL, purple_grap_cell_BL], name='purple_BL')

green_fuel_cell_BL = openmc.Cell(fill=f1, region=(    (inside_channel_1) ))
green_grap_cell_BL = openmc.Cell(fill=g1,region=(~inside_channel_1 | inside_channel_3 | inside_channel_2 | inside_channel_4))
green_BL = openmc.Universe(cells=[green_fuel_cell_BL, green_grap_cell_BL], name='green_BL')

yellow_fuel_cell_BL = openmc.Cell(fill=f1, region=(    (inside_channel_2) ))
yellow_grap_cell_BL = openmc.Cell(fill=g1,region=(~inside_channel_2 | inside_channel_1 | inside_channel_4 | inside_channel_3))
yellow_BL = openmc.Universe(cells=[yellow_fuel_cell_BL, yellow_grap_cell_BL], name='yellow_BL')

brown_fuel_cell_BL = openmc.Cell(fill=f1, region=(    (inside_channel_1 | inside_channel_3 | inside_channel_2) ))
brown_grap_cell_BL = openmc.Cell(fill=g1,region=(~inside_channel_1 & ~inside_channel_3 & ~inside_channel_2 | inside_channel_4))
brown_BL = openmc.Universe(cells=[brown_fuel_cell_BL, brown_grap_cell_BL], name='brown_BL')



### TOP LEFT QUARTER OF CORE
blue_fuel_cell_TL = openmc.Cell(fill=f1, region=(    (inside_channel_2 | inside_channel_3 | inside_channel_4) ))
blue_grap_cell_TL = openmc.Cell(fill=g1,region=(~inside_channel_2 & ~inside_channel_3 & ~inside_channel_4 | inside_channel_1))
blue_TL = openmc.Universe(cells=[blue_fuel_cell_TL, blue_grap_cell_TL], name='blue_TL')

purple_fuel_cell_TL = openmc.Cell(fill=f1, region=(    (inside_channel_3 | inside_channel_2) ))
purple_grap_cell_TL = openmc.Cell(fill=g1,region=(  ~inside_channel_3 & ~inside_channel_2  | inside_channel_1 | inside_channel_4  ))
purple_TL = openmc.Universe(cells=[purple_fuel_cell_TL, purple_grap_cell_TL], name='purple_TL')

green_fuel_cell_TL = openmc.Cell(fill=f1, region=(    (inside_channel_3) ))
green_grap_cell_TL = openmc.Cell(fill=g1,region=(~inside_channel_3 | inside_channel_3 | inside_channel_2 | inside_channel_1))
green_TL = openmc.Universe(cells=[green_fuel_cell_TL, green_grap_cell_TL], name='green_TL')

yellow_fuel_cell_TL = openmc.Cell(fill=f1, region=(    (inside_channel_2) ))
yellow_grap_cell_TL = openmc.Cell(fill=g1,region=(~inside_channel_2 | inside_channel_1 | inside_channel_4 | inside_channel_3))
yellow_TL = openmc.Universe(cells=[yellow_fuel_cell_TL, yellow_grap_cell_TL], name='yellow_TL')

brown_fuel_cell_TL = openmc.Cell(fill=f1, region=(    (inside_channel_1 | inside_channel_3 | inside_channel_2) ))
brown_grap_cell_TL = openmc.Cell(fill=g1,region=(~inside_channel_1 & ~inside_channel_3 & ~inside_channel_2 | inside_channel_4))
brown_TL = openmc.Universe(cells=[brown_fuel_cell_TL, brown_grap_cell_TL], name='brown_TL')


############################################
### FUEL ONLY CELL
############################################
fuelOnlyCell = openmc.Cell(fill=f1)
F = openmc.Universe()
F.add_cell(fuelOnlyCell)


############################################
### CONTROL ROD UNIT CELL
############################################
#### control rod 2.0 at rest in the top of the core
insertion_depth_minimum = 2.75 * 2.54 # 833ft,5.75in minus 833ft,3in ---> 2.75 inches ---> https://www.osti.gov/servlets/purl/1617123 page 15
upper_graphite_stringer_height = 67*2.54/2
lower_graphite_stringer_height = -67*2.54/2

# radial values
thimble_od = 5.08 # 5.08 outer diameter
thimble_thickness = 0.1651
poison_od = 1.08*2.54 # from ORNL 4123 page 37
poison_id = 0.840*2.54 # ORNL 4123  page 37
shell_thickness = 0.02*2.54 # inconel shell thickness that goes around each contrl element
ss_id = 5.0/8.0 * 2.54 # inner diameter of stainless steel threading
ss_od = 0.75 * 2.54  # outer diameter of stainless steel threading
cable_od_D = 0.3175

# axial values
poison_height = 1.325*2.54 # ORNL DWG 64-8816
can_height = 1.562 * 2.54   # ORNL DWG 64-8816
num_cans_high = 38 # number of stacks of cans
can_t = (can_height - poison_height) / 2

# axial planes for unit cell
poison_lower_plane = openmc.ZPlane(z0=can_t - can_height/2)
poison_upper_plane = openmc.ZPlane(z0=can_t+poison_height - can_height/2)
can_lower = openmc.ZPlane(z0=0.0 - can_height/2)
can_upper = openmc.ZPlane(z0=can_height - can_height/2)

# radial clinders for unit cell
# rod outer shell
thimble_outer = openmc.ZCylinder(r=thimble_od/2)
thimble_inner = openmc.ZCylinder(r=thimble_od/2 - thimble_thickness)
poison_outershell_od = openmc.ZCylinder(r=poison_od/2)
poison_outershell_id = openmc.ZCylinder(r=poison_od/2 - shell_thickness)
poison_innershell_od = openmc.ZCylinder(r=poison_id/2 + shell_thickness)
poison_innershell_id = openmc.ZCylinder(r=poison_id/2)
stainless_id = openmc.ZCylinder(r=ss_id/2)
stainless_od = openmc.ZCylinder(r=ss_od/2)
cable_od = openmc.ZCylinder(r=(cable_od_D/2) * 1.0)

# Cells for unit cell of 1 element
thimble_cell = openmc.Cell(region=(-thimble_outer & +thimble_inner & +can_lower & -can_upper), fill=inor)
air_cell_one = openmc.Cell(region = (-thimble_inner & +poison_outershell_od & +can_lower & -can_upper), fill=cellgas)
inc_shell_outer = openmc.Cell(region=(-poison_outershell_od & +poison_outershell_id & +can_lower & -can_upper), fill=inconel)
inc_shell_inner = openmc.Cell(region=(-poison_innershell_od & +poison_innershell_id & +can_lower & -can_upper), fill=inconel)
air_cell_two = openmc.Cell(region=(-poison_innershell_id & +stainless_od & +can_lower & -can_upper), fill=cellgas)
ss_filling = openmc.Cell(region=(-stainless_od & +stainless_id & +can_lower & -can_upper), fill=ss316)
air_cell_three = openmc.Cell(region=(-stainless_id & +cable_od & +can_lower & -can_upper), fill=cellgas)
cable_cell = openmc.Cell(region=(-cable_od & +can_lower & -can_upper), fill=inconel)


poison_shell = openmc.Cell(region=(-poison_outershell_id & +poison_innershell_od & +poison_lower_plane & -poison_upper_plane), fill=poison)
can_shell_upper = openmc.Cell(region=(-poison_outershell_id & +poison_innershell_od & +poison_upper_plane & -can_upper), fill=inconel)
can_shell_lower = openmc.Cell(region=(-poison_outershell_id & +poison_innershell_od & -poison_lower_plane & +can_lower), fill=inconel)

can_inf_fill = openmc.Cell(region=(+can_upper | - can_lower), fill=f1)
all_flibe = openmc.Cell(region=(+thimble_outer & +can_lower & -can_upper), fill=f1, name='all_flibe')
rodded_flibe_1 = openmc.Cell(region=(+thimble_outer & +can_lower & -can_upper), fill=f2, name='rodded_flibe_1')
rodded_flibe_2 = openmc.Cell(region=(+thimble_outer & +can_lower & -can_upper), fill=f2, name='rodded_flibe_2')

# univerese with only rod (no surrounding flibe)
R_uni = openmc.Universe(cells=[thimble_cell, air_cell_one, inc_shell_outer, inc_shell_inner, air_cell_two, ss_filling, air_cell_three, cable_cell,
                              poison_shell, can_shell_upper, can_shell_lower,
                              can_inf_fill, all_flibe])
# ,all_flibe

R_only_cell_1 = openmc.Cell(region=(-thimble_outer & +can_lower & -can_upper), fill=R_uni, name='R_only_cell_1')
R_only_cell_2 = openmc.Cell(region=(-thimble_outer & +can_lower & -can_upper), fill=R_uni, name='R_only_cell_2')
Rcell_1 = openmc.Universe(cells=[R_only_cell_1, rodded_flibe_1])
Rcell_2 = openmc.Universe(cells=[R_only_cell_2, rodded_flibe_2])

# Rcell_1.plot(pixels=(600,600), basis='xz',width=(8,8),origin=(0,0,0),color_by='cell')


############################################
### ROD LATTICE 1
############################################

# Section 2.1.25 from MSRE benchmark
# At the time of criticality, the position of the regulating rod at room temperature was 118.364 cm whereas
# the other two rods were 129.54 cm -> two rods=51 inch position and regulating rod = 4.4 in lower than passive rods. (4.4 inch insertion)

### ROD LATTICE NUMBER 1 FOR PLACEMENT IN CORE AND UPPER PLENUM
graphite_top_z = 67*2.54/2
graphite_bottom_z = -67*2.54/2
rodLowestLimit_r1 = graphite_bottom_z + 0.0### 3 ### originally 3 inhes but simply to 0.0 inches so bottom of rod is at bottom of the core.
rod_depth_r1 = rod_1_user_input_depth_inches*2.54 + 0.0*2.54#### + 2.75*2.54# (cm measured from top of graphite) aka insertion depth
lower_left_rod_z_r1 = graphite_top_z - rod_depth_r1

# lower rod planes and cells - part that goes into core
rodLatticeUpperZPlane_r1 = openmc.ZPlane(z0=graphite_top_z) # upper graphite
rodLatticeLowerZPlane_r1 = openmc.ZPlane(z0=lower_left_rod_z_r1) # minimum value of rod distance poison segment region
lower_rod_cap_plane_r1 = openmc.ZPlane(z0=lower_left_rod_z_r1-thimble_thickness) # plane for upper crod cap
lowestRodZPlane_r1 = openmc.ZPlane(z0=rodLowestLimit_r1) # place at bottom of core slightly above graphtie
lowerRodZPlanePlusThickness_r1 = openmc.ZPlane(rodLowestLimit_r1+thimble_thickness)
lowerCorePlaneForRod_r1 = openmc.ZPlane(z0=graphite_bottom_z)

lower_rod_shell_r1 = openmc.Cell(region=(-thimble_outer & +thimble_inner & +lowestRodZPlane_r1 & -rodLatticeLowerZPlane_r1), fill=inor_in_core, name='lower_rod_shell_r1') # shell of rod
lower_rod_cap_r1 = openmc.Cell(region=(-thimble_inner & +lower_rod_cap_plane_r1 & -rodLatticeLowerZPlane_r1), fill=inor_in_core, name='lower_rod_cap_r1') # rod cap upper part for air part of rod
lower_rod_air_r1 = openmc.Cell(region=(-thimble_inner & -lower_rod_cap_plane_r1 & +lowerRodZPlanePlusThickness_r1), fill=cellgas_in_core, name='lower_rod_air_r1') # air fill for lower empty part of rod
lower_rod_cap_lower_r1 = openmc.Cell(region=(-thimble_inner & +lowestRodZPlane_r1 & -lowerRodZPlanePlusThickness_r1), fill=inor_in_core, name='lower_rod_cap_lower_r1') # lower cap for the lower empty part of the rod
lower_rod_fuel_r1 = openmc.Cell(region=(+thimble_outer & +lowestRodZPlane_r1 & -rodLatticeLowerZPlane_r1), fill=f1, name='lower_rod_fuel_r1')  # fuel around rod
lower_unrodded_fuel_r1 = openmc.Cell(region=(+lowerCorePlaneForRod_r1 & -lowestRodZPlane_r1), fill=f1, name='lower_unrodded_fuel_r1') # fuel below rod - unrodded region

# upper rod planes and cells - part that goes above core in plenum - need to make continuous with part that goes in core
rod_plenum_upper_plane_r1 = openmc.ZPlane(z0=lower_left_rod_z_r1+num_cans_high*can_height)

rod_plenum_upper_inco_shell_r1 = openmc.Cell(region=(-thimble_outer & +thimble_inner & +rod_plenum_upper_plane_r1), fill=inor, name='rod_plenum_upper_inco_shell_r1')
rod_plenum_upper_outer_air_r1 = openmc.Cell(region=(+stainless_od & -thimble_inner & +rod_plenum_upper_plane_r1), fill=cellgas, name='rod_plenum_upper_outer_air_r1')
rod_plenum_upper_ss_wire_r1 = openmc.Cell(region=(-stainless_od & +stainless_id & +rod_plenum_upper_plane_r1), fill=ss316, name='rod_plenum_upper_ss_wire_r1')
rod_plenum_upper_inner_air_r1 = openmc.Cell(region=(-stainless_id & +cable_od & +rod_plenum_upper_plane_r1), fill=cellgas, name='rod_plenum_upper_inner_air_r1')
rod_plenum_upper_inner_cable_r1 = openmc.Cell(region=(-cable_od & +rod_plenum_upper_plane_r1), fill=inor, name='rod_plenum_upper_inner_cable_r1')
rod_plenum_upper_fuel_r1 = openmc.Cell(region=(+thimble_outer & +rod_plenum_upper_plane_r1), fill=f1, name='rod_plenum_upper_fuel_r1')

rodLattice_r1=openmc.RectLattice()
rodLattice_r1.lower_left = (-L/2,-L/2,lower_left_rod_z_r1)
rodLattice_r1.pitch = (L,L,can_height)
rodLattice_r1_fuel_cell = openmc.Cell(fill=f1, name='rodLattice_r1_fuel_cell')
rodLattice_r1.outer = openmc.Universe(cells=[rodLattice_r1_fuel_cell,])
rodLattice_r1.universes = [[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],
                           [[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],
                           [[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],
                           [[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]],[[Rcell_1]]]




lowerRodLatticeCell_r1 = openmc.Cell(region=(+rodLatticeLowerZPlane_r1 & -rodLatticeUpperZPlane_r1), fill=rodLattice_r1)
upperRodLatticeCell_r1 = openmc.Cell(region=(-rod_plenum_upper_plane_r1 & +rodLatticeUpperZPlane_r1), fill=rodLattice_r1)



RLT1 = openmc.Universe(cells=[lowerRodLatticeCell_r1, upperRodLatticeCell_r1,
                              lower_rod_shell_r1, lower_rod_cap_r1, lower_rod_air_r1, lower_rod_cap_lower_r1,
                              lower_rod_fuel_r1, lower_unrodded_fuel_r1,
                              rod_plenum_upper_inco_shell_r1, rod_plenum_upper_outer_air_r1, rod_plenum_upper_ss_wire_r1, rod_plenum_upper_inner_air_r1, rod_plenum_upper_inner_cable_r1,
                              rod_plenum_upper_fuel_r1])



# RLT1.plot(pixels=(2000,2000), basis='xz',width=(8,500),origin=(0,0,0),color_by='material', colors={cellgas:'yellow', inconel:'black', f1:'blue', poison:'pink', ss316:"green", g1:'teal'})
# rodLatticeUni.plot(pixels=(300,300), basis='xy',width=(6,6),origin=(0,0,can_height/2),color_by='material', colors={cellgas:'yellow', inconel:'black', f1:'blue', poison:'pink', ss316:"green", g1:'teal'})




############################################
### ROD LATTICE 2
############################################
### ROD LATTICE NUMBER 2 FOR PLACEMENT IN CORE AND UPPER PLENUM
graphite_top_z = 67*2.54/2
graphite_bottom_z = -67*2.54/2
rodLowestLimit_r2 = graphite_bottom_z + 0.0 ### 3 ## should be 3 but set to 0 to simplify
rod_depth_r2 = rod_2_user_input_depth_inches*2.54 + 0.0*2.54 #### 2.75*2.54 # (cm measured from top of graphite) aka insertion depth
lower_left_rod_z_r2 = graphite_top_z - rod_depth_r2

# lower rod planes and cells - part that goes into core
rodLatticeUpperZPlane_r2 = openmc.ZPlane(z0=graphite_top_z) # upper graphite
rodLatticeLowerZPlane_r2 = openmc.ZPlane(z0=lower_left_rod_z_r2) # minimum value of rod distance poison segment region
lower_rod_cap_plane_r2 = openmc.ZPlane(z0=lower_left_rod_z_r2-thimble_thickness) # plane for upper crod cap
lowestRodZPlane_r2 = openmc.ZPlane(z0=rodLowestLimit_r2) # place at bottom of core slightly above graphtie
lowerRodZPlanePlusThickness_r2 = openmc.ZPlane(rodLowestLimit_r2+thimble_thickness)
lowerCorePlaneForRod_r2 = openmc.ZPlane(z0=graphite_bottom_z)

lower_rod_shell_r2 = openmc.Cell(region=(-thimble_outer & +thimble_inner & +lowestRodZPlane_r2 & -rodLatticeLowerZPlane_r2), fill=inor_in_core, name='lower_rod_shell_r2') # shell of rod
lower_rod_cap_r2 = openmc.Cell(region=(-thimble_inner & +lower_rod_cap_plane_r2 & -rodLatticeLowerZPlane_r2), fill=inor_in_core, name='lower_rod_cap_r2') # rod cap upper part for air part of rod
lower_rod_air_r2 = openmc.Cell(region=(-thimble_inner & -lower_rod_cap_plane_r2 & +lowerRodZPlanePlusThickness_r2), fill=cellgas_in_core, name='lower_rod_air_r2') # air fill for lower empty part of rod
lower_rod_cap_lower_r2 = openmc.Cell(region=(-thimble_inner & +lowestRodZPlane_r2 & -lowerRodZPlanePlusThickness_r2), fill=inor_in_core, name='lower_rod_cap_lower_r2') # lower cap for the lower empty part of the rod
lower_rod_fuel_r2 = openmc.Cell(region=(+thimble_outer & +lowestRodZPlane_r2 & -rodLatticeLowerZPlane_r2), fill=f1, name='lower_rod_fuel_r2')  # fuel around rod
lower_unrodded_fuel_r2 = openmc.Cell(region=(+lowerCorePlaneForRod_r2 & -lowestRodZPlane_r2), fill=f1, name='lower_unrodded_fuel_r2') # fuel below rod - unrodded region

# upper rod planes and cells - part that goes above core in plenum - need to make continuous with part that goes in core
rod_plenum_upper_plane_r2 = openmc.ZPlane(z0=lower_left_rod_z_r2+num_cans_high*can_height)

rod_plenum_upper_inco_shell_r2 = openmc.Cell(region=(-thimble_outer & +thimble_inner & +rod_plenum_upper_plane_r2), fill=inor, name='rod_plenum_upper_inco_shell_r2')
rod_plenum_upper_outer_air_r2 = openmc.Cell(region=(+stainless_od & -thimble_inner & +rod_plenum_upper_plane_r2), fill=cellgas, name='rod_plenum_upper_outer_air_r2')
rod_plenum_upper_ss_wire_r2 = openmc.Cell(region=(-stainless_od & +stainless_id & +rod_plenum_upper_plane_r2), fill=ss316, name='rod_plenum_upper_ss_wire_r2')
rod_plenum_upper_inner_air_r2 = openmc.Cell(region=(-stainless_id & +cable_od & +rod_plenum_upper_plane_r2), fill=cellgas, name='rod_plenum_upper_inner_air_r2')
rod_plenum_upper_inner_cable_r2 = openmc.Cell(region=(-cable_od & +rod_plenum_upper_plane_r2), fill=inconel, name='rod_plenum_upper_inner_cable_r2')
rod_plenum_upper_fuel_r2 = openmc.Cell(region=(+thimble_outer & +rod_plenum_upper_plane_r2), fill=f1, name='rod_plenum_upper_fuel_r2')

rodLattice_r2=openmc.RectLattice()
rodLattice_r2.lower_left = (-L/2,-L/2,lower_left_rod_z_r2)
rodLattice_r2.pitch = (L,L,can_height)
rodLattice_r2_fuel_cell = openmc.Cell(fill=f1, name='rodLattice_r2_fuel_cell')
rodLattice_r2.outer = openmc.Universe(cells=[rodLattice_r2_fuel_cell,])
rodLattice_r2.universes = [[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],
                           [[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],
                           [[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],
                           [[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]],[[Rcell_2]]]




lowerRodLatticeCell_r2 = openmc.Cell(region=(+rodLatticeLowerZPlane_r2 & -rodLatticeUpperZPlane_r2), fill=rodLattice_r2)
upperRodLatticeCell_r2 = openmc.Cell(region=(-rod_plenum_upper_plane_r2 & +rodLatticeUpperZPlane_r2), fill=rodLattice_r2)



RLT2 = openmc.Universe(cells=[lowerRodLatticeCell_r2, upperRodLatticeCell_r2,
                              lower_rod_shell_r2, lower_rod_cap_r2, lower_rod_air_r2, lower_rod_cap_lower_r2,
                              lower_rod_fuel_r2, lower_unrodded_fuel_r2,
                              rod_plenum_upper_inco_shell_r2, rod_plenum_upper_outer_air_r2, rod_plenum_upper_ss_wire_r2, rod_plenum_upper_inner_air_r2, rod_plenum_upper_inner_cable_r2,
                              rod_plenum_upper_fuel_r2])



# RLT2.plot(pixels=(2000,2000), basis='xz',width=(8,500),origin=(0,0,0),color_by='material', colors={cellgas:'yellow', inconel:'black', f1:'blue', poison:'pink', ss316:"green", g1:'teal'})



############################################
### BASKET MODELLING
############################################
##### BASKET MODELLING

inor_od = 5.08 # 5.08 outer diameter
inor_t = 0.079
inor_id = inor_od - inor_t*2

inor_id_surface = openmc.ZCylinder(r=inor_id/2.0)
inor_od_surface = openmc.ZCylinder(r=inor_od/2.0)


inor_cell_basket = openmc.Cell(region=(-inor_od_surface & +inor_id_surface),fill=inor)

gWidth = 0.6 * 2.54
gHeight = 0.5 * 2.54

gLower = openmc.YPlane(y0=-gHeight/2)
gUpper = openmc.YPlane(y0=gHeight/2)
gLeft = openmc.XPlane(x0=-gWidth)
gRight = openmc.XPlane(x0=gWidth)

gUpper2 = openmc.YPlane(y0=gHeight/2 + gHeight)
gLeft2 = openmc.XPlane(x0=-gWidth/2)
gRight2 = openmc.XPlane(x0=gWidth/2)

flibe_in_basket = openmc.Cell(region=(-inor_id_surface), fill=f1)
graphite_block_1 = openmc.Cell(region=(+gLower & -gUpper & +gLeft & -gRight) , fill=g1)
graphite_block_2 = openmc.Cell(region=(+gUpper & +gLeft2 & -gRight2 & -gUpper2), fill=g1)

inor_cyl_rad = 0.25 / 2 * 2.54 # quarter nich in diameter
inor_cylinder_1 = openmc.ZCylinder(r=inor_cyl_rad, x0=-gWidth-inor_cyl_rad, y0=0)
inor_cylinder_2 = openmc.ZCylinder(r=inor_cyl_rad, x0=-gWidth/2-inor_cyl_rad, y0=gHeight)
inor_cylinder_3 = openmc.ZCylinder(r=inor_cyl_rad, x0= gWidth/2+inor_cyl_rad, y0=gHeight)
inor_cylinder_4 = openmc.ZCylinder(r=inor_cyl_rad, x0=+gWidth+inor_cyl_rad, y0=0)
inor_cylinder_5 = openmc.ZCylinder(r=inor_cyl_rad, x0=+gWidth/2, y0=-inor_cyl_rad-gHeight/2)
inor_cylinder_6 = openmc.ZCylinder(r=inor_cyl_rad, x0=-gWidth/2, y0=-inor_cyl_rad-gHeight/2)


inor_cyl_1 = openmc.Cell(region=(-inor_cylinder_1 & -gLeft), fill=inor)
inor_cyl_2 = openmc.Cell(region=(-inor_cylinder_2 & -gLeft2), fill=inor)
inor_cyl_3 = openmc.Cell(region=(-inor_cylinder_3 & +gRight2), fill=inor)
inor_cyl_4 = openmc.Cell(region=(-inor_cylinder_4 & +gRight), fill=inor)
inor_cyl_5 = openmc.Cell(region=(-inor_cylinder_5 & -gLower), fill=inor)
inor_cyl_6 = openmc.Cell(region=(-inor_cylinder_6 & -gLower), fill=inor)

region_flibe = ~inor_cyl_1.region & ~inor_cyl_2.region & ~inor_cyl_3.region & ~inor_cyl_4.region & ~inor_cyl_5.region & ~inor_cyl_6.region & ~graphite_block_1.region & ~graphite_block_2.region & -inor_id_surface
flibe_in_basket = openmc.Cell(region=region_flibe, fill=f1)

all_flibe_around_basket = openmc.Cell(region=(+inor_id_surface), fill=f1)

K = openmc.Universe(cells=[inor_cell_basket, all_flibe_around_basket, graphite_block_1,graphite_block_2, inor_cyl_1, inor_cyl_4, inor_cyl_2, inor_cyl_3, inor_cyl_5, inor_cyl_6, flibe_in_basket])
# K.plot(width=(5,5), color_by='cell')


############################################
### SOFT GRAPHITE AND FUEL CHANNELS THAT GO AROUND CONTROL RODS
############################################
########### universe for the soft graphite channels around the control rods that are placed in the core
# inside_channel_1 = ( (-c1 | -c2) & (+p6 & -p2) ) | (+p9 & -p10 & +p6 & -p2)
# inside_channel_2 = ( (-c3 | -c4) & (+p7 & -p3) ) | (+p12 & -p11 & +p7 & -p3)
# inside_channel_3 = ( (-c5 | -c6) & (+p4 & -p8) ) | (+p9 & -p10 & +p4 & -p8)
# inside_channel_4 = ( (-c7 | -c8) & (+p1 & -p5) ) | (+p12 & -p11 & +p1 & -p5)
rad_channel_special = 3.1992
#### channels placed vertically (cylinders on left and right sides)

left_shifted_cyl = openmc.ZCylinder(r=rad_channel_special, x0=-L)
right_shifted_cyl = openmc.ZCylinder(r=rad_channel_special,x0=L)

vert_g_cell = openmc.Cell(region=(~inside_channel_1 & ~inside_channel_3 & +left_shifted_cyl & + right_shifted_cyl),fill=g1)
vert_flibe_cells = openmc.Cell(region=(-left_shifted_cyl | -right_shifted_cyl | inside_channel_1 | inside_channel_3),fill=f1)

V = openmc.Universe(cells=[vert_g_cell,vert_flibe_cells,])

# V.plot(width=(15,15))


#### channels placed horizontally
up_shifted_cyl = openmc.ZCylinder(r=rad_channel_special, y0=L)
down_shifted_cyl = openmc.ZCylinder(r=rad_channel_special,y0=-L)
hori_g_cell = openmc.Cell(region=(~inside_channel_2 & ~inside_channel_4 & +up_shifted_cyl & +down_shifted_cyl),fill=g1)
hori_flibe_cells = openmc.Cell(region=(inside_channel_2 | inside_channel_4 | -up_shifted_cyl | -down_shifted_cyl),fill=f1)

H = openmc.Universe(cells=[hori_flibe_cells, hori_g_cell])


#### channels placed horizontally on TOP of the thing
hori_g_cell_top = openmc.Cell(region=(~inside_channel_2 & ~inside_channel_4 & ~inside_channel_1 & +down_shifted_cyl),fill=g1)
hori_flibe_cell_top = openmc.Cell(region=(inside_channel_2 | inside_channel_4 | inside_channel_1 | -down_shifted_cyl),fill=f1)
T = openmc.Universe(cells=[hori_g_cell_top, hori_flibe_cell_top])

#### channels placed horizontally on BOTTOM of the thing
hori_g_cell_bot = openmc.Cell(region=(~inside_channel_2 & ~inside_channel_4 & ~inside_channel_3 & +up_shifted_cyl),fill=g1)
hori_flibe_cell_bot = openmc.Cell(region=(inside_channel_2 | inside_channel_4 | inside_channel_3 | -up_shifted_cyl),fill=f1)
B = openmc.Universe(cells=[hori_g_cell_bot, hori_flibe_cell_bot])


#### channels placed vertically on LEFT of the thing
vert_g_cell_left = openmc.Cell(region=(~inside_channel_1 & ~inside_channel_3 & ~inside_channel_4 & + right_shifted_cyl),fill=g1)
vert_flibe_cell_left = openmc.Cell(region=(-right_shifted_cyl | inside_channel_1 | inside_channel_3 | inside_channel_4),fill=f1)
M = openmc.Universe(cells=[vert_flibe_cell_left, vert_g_cell_left])

#### channels placed vertically on RIGHT of the thing
vert_g_cell_right = openmc.Cell(region=(~inside_channel_1 & ~inside_channel_3 & ~inside_channel_2 & + left_shifted_cyl),fill=g1)
vert_flibe_cell_right = openmc.Cell(region=(-left_shifted_cyl | inside_channel_1 | inside_channel_3 | inside_channel_2),fill=f1)
N = openmc.Universe(cells=[vert_flibe_cell_right, vert_g_cell_right])

# R.plot(width=(10,10))



############################################
### GEOMETRY ALL DUMP
############################################

lattice=openmc.RectLattice()
lattice.lower_left = (-L*15.5,-L*15.5)
lattice.pitch = (L,L)
lattice.outer = G
lattice.universes = [
                     [G,G,  G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,  G,G],
                     [G,G,  G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,  G,G],

                     [G,G,  G,G,G,G,G,G,G,G,G,green_TL,purple_TL,blue_TL,u  ,u  ,u,blue_TR,purple_TR,green_TR,G,G,G,G,G,G,G,G,G,  G,G], # 13 | 1
                     [G,G,  G,G,G,G,G,G,G,purple_TL,u,u,u,u,u  ,u  ,u,u,u,u,u,purple_TR,G,G,G,G,G,G,G,  G,G], # 12| 2
                     [G,G,  G,G,G,G,G,purple_TL,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,purple_TR,G,G,G,G,G,  G,G], # 11| 3
                     [G,G,  G,G,G,G,purple_TL,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,purple_TR,G,G,G,G,  G,G], # 10| 4
                     [G,G,  G,G,G,purple_TL,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,purple_TR,G,G,G,  G,G], # 9| 5
                     [G,G,  G,G,purple_TL,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,purple_TR,G,G,  G,G], # 8| 6
                     [G,G,  G,G,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,G,G,  G,G], # 7| 7
                     [G,G,  G,purple_TL,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,purple_TR,G,  G,G], # 6| 8
                     [G,G,  G,u,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,u,G,  G,G], # 5| 9
                     [G,G,  yellow_TL,u,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,u,yellow_TR,  G,G], # 4 | 10
                     [G,G,  purple_TL,u,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,u,purple_TR,  G,G], # 3| 11
                     [G,G,  brown_TL,u,u,u,u,u,u,u,u,u,u,u,T  ,u  ,T,u,u,u,u,u,u,u,u,u,u,u,brown_TR,  G,G], # 2| 12
                     [G,G,  u,u,u,u,u,u,u,u,u,u,u,M,RLT2  ,V  ,RLT1,N,u,u,u,u,u,u,u,u,u,u,u,  G,G], # 1| 13

                     [G,G,  u,u,u,u,u,u,u,u,u,u,u,u,H  ,u  ,B,u,u,u,u,u,u,u,u,u,u,u,u,  G,G],         # 14

                     [G,G,  u,u,u,u,u,u,u,u,u,u,u,M,RLT1  ,N     ,u,u,u,u,u,u,u,u,u,u,u,u,u,  G,G], # 1 --- K s the basket
                     [G,G,  brown_BL,u,u,u,u,u,u,u,u,u,u,u,B  ,u  ,u,u,u,u,u,u,u,u,u,u,u,u,brown_BR,  G,G], # 2
                     [G,G,  purple_BL,u,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,u,purple_BR,  G,G], # 3
                     [G,G,  yellow_BL,u,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,u,yellow_BR,  G,G], # 4
                     [G,G,  G,u,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,u,G,  G,G], # 5
                     [G,G,  G,purple_BL,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,purple_BR,G,  G,G], # 6
                     [G,G,  G,G,u,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,u,G,G,  G,G], # 7
                     [G,G,  G,G,purple_BL,u,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,u,purple_BR,G,G,  G,G], # 8
                     [G,G,  G,G,G,purple_BL,u,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,u,purple_BR,G,G,G,  G,G], # 9
                     [G,G,  G,G,G,G,purple_BL,u,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,u,purple_BR,G,G,G,G,  G,G], # 10
                     [G,G,  G,G,G,G,G,purple_BL,u,u,u,u,u,u,u  ,u  ,u,u,u,u,u,u,u,purple_BR,G,G,G,G,G,  G,G], # 11
                     [G,G,  G,G,G,G,G,G,G,purple_BL,u,u,u,u,u  ,u  ,u,u,u,u,u,purple_BR,G,G,G,G,G,G,G,  G,G], # 12
                     [G,G,  G,G,G,G,G,G,G,G,G,green_BL,purple_BL,blue_BL,u  ,u  ,u,blue_BR,purple_BR,green_BR,G,G,G,G,G,G,G,G,G,  G,G], # 13

                     [G,G,  G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,  G,G],
                     [G,G,  G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,  G,G]
]

# print(lattice)

# core surfaces for graphite regions:
core_upper_z = 67*2.54/2
core_lower_z = -67*2.54/2
upper_core_surface = openmc.ZPlane(z0=core_upper_z)
lower_core_surface = openmc.ZPlane(z0=core_lower_z)
core_height_bounds = (-upper_core_surface & +lower_core_surface)


# surface for cylinders around the core
outer_core_surface = openmc.ZCylinder(r=70.285) ### TODO needs correction ?
inner_can_surface =  openmc.ZCylinder(r=71.097)
outer_can_surface =  openmc.ZCylinder(r=71.737)
inner_vessel_surface = openmc.ZCylinder(r=74.299)
outer_vessel_surface = openmc.ZCylinder(r=76.862)   # , boundary_type = 'vacuum')
core = openmc.Cell(fill=lattice, region=-outer_core_surface & (-upper_core_surface & +lower_core_surface))

# TODO FIX FILLS HERE BY ADDING CORRECT MATERIALS
upcomer_cell = openmc.Cell(region=(+outer_core_surface & -inner_can_surface & (-upper_core_surface & +lower_core_surface)), fill=f2)
can_cell = openmc.Cell(region=(-outer_can_surface & +inner_can_surface & (-upper_core_surface & +lower_core_surface)), fill=inor)
downcomer_cell = openmc.Cell(region=(-inner_vessel_surface & +outer_can_surface & (-upper_core_surface & +lower_core_surface)), fill=f2)
vessel_cell = openmc.Cell(region=(-outer_vessel_surface & +inner_vessel_surface & (-upper_core_surface & +lower_core_surface)), fill=inor)


# universes for area surrounding the main core
upcomer_u = openmc.Universe(cells=[upcomer_cell])
downcomer_u = openmc.Universe(cells=[downcomer_cell])
vessel_u = openmc.Universe(cells=[vessel_cell])
can_u = openmc.Universe(cells=[can_cell])


# 3 control rods in upper part of the core:
upperRightCylinder = openmc.ZCylinder(r=thimble_od/2, x0=L, y0=L)
upperLeftCylinder = openmc.ZCylinder(r=thimble_od/2, x0=-L, y0=L)
lowerLeftCylinder = openmc.ZCylinder(r=thimble_od/2, x0=-L, y0=-L)
# squares around control rods in plena (need to make crods square so that they match the griffin meshing)
xFarLeft = openmc.XPlane(x0 = -1.5*L)
xFarRight = openmc.XPlane(x0 = 1.5*L)
xCloseLeft = openmc.XPlane(x0 = -0.5*L)
xCloseRight = openmc.XPlane(x0 = 0.5*L)
yFarTop = openmc.YPlane(y0 = 1.5*L)
yCloseTop = openmc.YPlane(y0 = 0.5*L)
yFarBottom = openmc.YPlane(y0 = -1.5*L)
yCloseBottom = openmc.YPlane(y0 = -0.5*L)

# regions for 3 different plena control rods.
upperLeftRegionFill = (+xFarLeft & -xCloseLeft & -yFarTop & +yCloseTop)
upperRightRegionFill = (-xFarRight & +xCloseRight & -yFarTop & +yCloseTop)
lowerLeftRegionFill = (+xFarLeft & -xCloseLeft & +yFarBottom & -yCloseBottom)

# upperRight = openmc.Cell(region=(-upperRightCylinder), fill=RLT1, name='upperRightRod')
upperRight = openmc.Cell(region=upperRightRegionFill, fill=RLT1, name='upperRightRod')
upperRight.translation=[L, L, 0]
# upperLeft = openmc.Cell(region=(-upperLeftCylinder), fill=RLT2, name='upperLeftRod')
upperLeft = openmc.Cell(region=(upperLeftRegionFill), fill=RLT2, name='upperLeftRod')
upperLeft.translation=[-1*L, L, 0]
# lowerLeft = openmc.Cell(region=(-lowerLeftCylinder), fill=RLT1, name='lowerLeftRod')
lowerLeft = openmc.Cell(region=(lowerLeftRegionFill), fill=RLT1, name='lowerLeftRod')
lowerLeft.translation=[-1*L, -1*L, 0]

# fuelAroundFourRodsRegion = ~(-lowerLeftCylinder | -upperRightCylinder | -upperLeftCylinder)
fuelAroundFourRodsRegion = ~(lowerLeftRegionFill | upperRightRegionFill | upperLeftRegionFill)

allButRods = openmc.Cell(region=fuelAroundFourRodsRegion, fill=f2, name="allButRods")

threeRods = openmc.Universe(cells=[upperLeft, upperRight, lowerLeft, allButRods])


### Upper plenum creation
plenum_height = 22.87953991 # based on elevations of MSRE (height of inner plena top shell minus top of graphite)
vessel_thickness = outer_vessel_surface.r - inner_vessel_surface.r # from vessel inner and outer surfaces ( approx 1 inch )
print("Vessel thickness is:",vessel_thickness)
print("Upper plena height is:", plenum_height)
r_outer = outer_vessel_surface.r # outside head radius -> simply the radius of the outside head. (radius of outer vessel wall)
r_crown = 2 * r_outer# inner crown radius  -> from sphere center to inner vessel top
r_knuckle = 0.03*r_outer # need to solve for this one # inside knuckle radius -> DETERMINED by solving:     (Rc - a)^2 = (Rc+t-h)^2 + (Router - t - a)^2

angleTheta = np.arccos(   (r_crown + vessel_thickness - plenum_height) / (r_crown - r_knuckle)  )
critical_radius_upper = (r_crown)*np.sin( angleTheta )
print('Upper Critical Radius = ', critical_radius_upper)

plenum_sphere_center = core_upper_z + plenum_height - r_crown

plenum_start = core_upper_z # use upper_core_surface for this surface
plenum_end = core_upper_z + plenum_height

upper_plena_vessel_outer_wall = openmc.Sphere(x0=0.0, y0=0.0, z0=plenum_sphere_center, r=r_crown) # , boundary_type = 'vacuum')
upper_plena_vessel_inner_wall = openmc.Sphere(x0=0.0, y0=0.0, z0=plenum_sphere_center, r=r_crown-vessel_thickness)
upper_vessel_crit_rad_cyl = openmc.ZCylinder(r=critical_radius_upper)

zTorus_inner = openmc.ZTorus(x0=0,y0=0, z0=core_upper_z, a=r_outer-r_knuckle-vessel_thickness, c=r_knuckle, b=r_knuckle)
zTorus_outer = openmc.ZTorus(x0=0,y0=0, z0=core_upper_z, a=r_outer-r_knuckle-vessel_thickness, c=r_knuckle+vessel_thickness, b=r_knuckle+vessel_thickness) # , boundary_type = 'vacuum')

upper_plenum_region = (-upper_vessel_crit_rad_cyl & +zTorus_inner & -zTorus_outer & -upper_plena_vessel_outer_wall & +upper_plena_vessel_inner_wall) | (+zTorus_inner & -zTorus_outer & +upper_core_surface & +upper_vessel_crit_rad_cyl & -outer_vessel_surface) | (-upper_plena_vessel_outer_wall & +upper_plena_vessel_inner_wall & +zTorus_outer & +upper_core_surface & -upper_vessel_crit_rad_cyl)
upper_plenum_fuel_region = (-zTorus_inner & +upper_core_surface) | (+upper_core_surface & -upper_plena_vessel_inner_wall & -upper_vessel_crit_rad_cyl)

upper_plena_shell_cell = openmc.Cell(region=upper_plenum_region, fill=inor)
upper_plenum_fuel_cell = openmc.Cell(region=upper_plenum_fuel_region, fill=threeRods)

### lower plenum creation
r_crown = 1.6*r_outer
r_knuckle = 0.308*r_outer
plenum_height = 39.76404
print("Lower plena height is:", plenum_height)
critical_radius_upper = (r_crown)*np.sin( np.arccos(   (r_crown + vessel_thickness - plenum_height) / (r_crown - r_knuckle)  ) )
upper_vessel_crit_rad_cyl = openmc.ZCylinder(r=critical_radius_upper)

print('Lower Critical Radius = ', critical_radius_upper)


lower_plena_vessel_outer_wall = openmc.Sphere(x0=0.0, y0=0.0, z0=core_lower_z - plenum_height + r_crown, r=r_crown) # , boundary_type = 'vacuum')
lower_plena_vessel_inner_wall = openmc.Sphere(x0=0.0, y0=0.0, z0=core_lower_z - plenum_height + r_crown, r=r_crown-vessel_thickness)
zTorus_lower_inner = openmc.ZTorus(x0=0,y0=0, z0=core_lower_z, a=r_outer-r_knuckle-vessel_thickness, c=r_knuckle, b=r_knuckle)
zTorus_lower_outer = openmc.ZTorus(x0=0,y0=0, z0=core_lower_z, a=r_outer-r_knuckle-vessel_thickness, c=r_knuckle+vessel_thickness, b=r_knuckle+vessel_thickness) # , boundary_type = 'vacuum')

lower_plenum_region = (-upper_vessel_crit_rad_cyl & +zTorus_lower_inner & -zTorus_lower_outer & -lower_plena_vessel_outer_wall & +lower_plena_vessel_inner_wall) | (+zTorus_lower_inner & -zTorus_lower_outer & -lower_core_surface & +upper_vessel_crit_rad_cyl & -outer_vessel_surface) | (-lower_plena_vessel_outer_wall & +lower_plena_vessel_inner_wall & +zTorus_lower_outer & -lower_core_surface & -upper_vessel_crit_rad_cyl)
lower_plenum_fuel_region = (-zTorus_lower_inner & -lower_core_surface) | (-lower_core_surface & -lower_plena_vessel_inner_wall & -upper_vessel_crit_rad_cyl)

lower_plena_fuel_cell = openmc.Cell(region=lower_plenum_fuel_region, fill=mixedLowerPlenaMat)

lower_plena_shell_cell = openmc.Cell(region=lower_plenum_region, fill=inor)


### THERMAL SHIELD
ts_height = 12.5 * 12 * 2.54
ts_ID = 7.5 * 2.54 * 12
ts_OD = 10.4 * 2.54 * 12
ts_thickness = ts_OD/2 - ts_ID/2
insulation_thickness = 6*2.54
ss_thickness = 1 * 2.54
axial_offset = ts_thickness/2


ts_od_cyl = openmc.ZCylinder(r=ts_OD/2)
ts_id_cyl = openmc.ZCylinder(r=ts_ID/2)
ts_ss_outer = openmc.ZCylinder(r=ts_OD/2 + ss_thickness, boundary_type = 'vacuum') # OK
ts_ss_inner = openmc.ZCylinder(r=ts_ID/2 - ss_thickness)
ts_ins_inner = openmc.ZCylinder(r=ts_ID/2 - ss_thickness - insulation_thickness)

ts_upper_plane = openmc.ZPlane(z0=ts_height/2                                                                       +axial_offset, boundary_type = 'vacuum', name='t_shield_upper_plane') # OK
ts_upper_ss_plate = openmc.ZPlane(z0=ts_height/2 + ss_thickness                                                     +axial_offset)
ts_lower_plane = openmc.ZPlane(z0=-ts_height/2                                                                      +axial_offset, boundary_type = 'vacuum', name='t_shield_lower_plane') # OK
ts_hori_lower_part = openmc.ZPlane(z0=ts_height/2-ts_thickness                                                      +axial_offset) # horizontal part - lower plane
ts_hori_lower_part_lower_plate = openmc.ZPlane(z0=ts_height/2-ts_thickness-ss_thickness                             +axial_offset) # lower plate bottom plane of the horizontal part
insulation_plane_upper = openmc.ZPlane(z0=ts_height/2-ts_thickness- ss_thickness - insulation_thickness             +axial_offset)
insulation_plane_lower = openmc.ZPlane(z0=-ts_height/2 + insulation_thickness             +axial_offset)

ts_cell = openmc.Cell(region = (-ts_od_cyl & +ts_id_cyl & -ts_upper_plane & +ts_lower_plane), fill=tshield)
ts_upper_stainless = openmc.Cell(region=(-ts_upper_ss_plate & +ts_upper_plane & -ts_od_cyl), fill=ss316)
ts_outer_stainless = openmc.Cell(region=(-ts_upper_ss_plate & -ts_ss_outer & +ts_od_cyl & +ts_lower_plane), fill=ss316)
ts_hori_cyl = openmc.Cell(region = (-ts_id_cyl & +ts_hori_lower_part & -ts_upper_plane), fill=tshield) # upper horizontal part of the thermal shield
ts_hori_cyl_lower_plate = openmc.Cell(region=(-ts_id_cyl & +ts_hori_lower_part_lower_plate & -ts_hori_lower_part), fill=ss316)
ts_inner_plate = openmc.Cell(region = (-ts_hori_lower_part_lower_plate & +ts_lower_plane & -ts_id_cyl & +ts_ss_inner), fill=ss316)


insulation_cell_upper = openmc.Cell(region = (-ts_hori_lower_part_lower_plate & +insulation_plane_upper & -ts_ss_inner), fill=insulation_mat)
insulation_inner_wall = openmc.Cell(region = (+ts_lower_plane & -insulation_plane_upper & -ts_ss_inner & +ts_ins_inner), fill=insulation_mat)
insulation_lower = openmc.Cell(region = (+ts_lower_plane & -insulation_plane_lower & -ts_ins_inner), fill=insulation_mat)

### external vessel fill gas
gas_region_1 = (-ts_ins_inner & +insulation_plane_lower & -insulation_plane_upper & +outer_vessel_surface)
gas_region_2 = (~upper_plenum_fuel_region & +upper_core_surface & ~upper_plenum_region & -insulation_plane_upper & -outer_vessel_surface)
gas_region_3 = (~lower_plenum_fuel_region & -lower_core_surface & ~lower_plenum_region & +insulation_plane_lower & -outer_vessel_surface)

gas_cell = openmc.Cell(region = (gas_region_1 | gas_region_2 | gas_region_3), fill=cellgas)

upcomer_cell.name = "upcomer_cell"
downcomer_cell.name = "downcomer_cell"
vessel_cell.name = "vessel_cell"
can_cell.name = "can_cell"
upper_plena_shell_cell.name = "upper_plena_shell_cell"
upper_plenum_fuel_cell.name = "upper_plena_fuel_cell"  ### UPPER PLENA WITH CONTROL RODS INCLUDED
lower_plena_shell_cell.name = "lower_plena_shell_cell" ### LOWER PLENA SHELL
lower_plena_fuel_cell.name = "lower_plena_fuel_cell"   ### LOWER PLENA FUEL
ts_cell.name = "ts_cell"
ts_upper_stainless.name = "ts_upper_stainless"
ts_outer_stainless.name = "ts_outer_stainless"
ts_hori_cyl_lower_plate.name = "ts_hori_cyl_lower_plate"
ts_inner_plate.name = "ts_inner_plate"
ts_hori_cyl.name = "ts_hori_cyl"
insulation_cell_upper.name = "insulation_cell_upper"
insulation_inner_wall.name = "insulation_inner_wall"
insulation_lower.name = "insulation_lower"
gas_cell.name = "gas_cell"

# crod related stuff / cells
rodLattice_r1_fuel_cell.name = 'rodLattice_r1_fuel_cell'
rodLattice_r2_fuel_cell.name = 'rodLattice_r2_fuel_cell'
lower_unrodded_fuel_r1.name = 'lower_unrodded_fuel_r1' # this is below the physical control rod in the small region without fuels
lower_unrodded_fuel_r2.name = 'lower_unrodded_fuel_r2' # same as text above
lower_rod_fuel_r1.name = 'lower_rod_fuel_r1'           # lower part of the crod with fuel in it - for part with fuel and inco and air only - rod 1
lower_rod_fuel_r2.name = 'lower_rod_fuel_r2'           # lower part of the crod with fuel in it - for part with fuel and inco and air only - rod 2
lower_rod_air_r1.name = 'lower_rod_air_r1'             # # lower part of the crod with AIR in it - for part with fuel and inco and air only - rod 1
lower_rod_air_r2.name = 'lower_rod_air_r2'             # # lower part of the crod with AIR in it - for part with fuel and inco and air only - rod 2
lower_rod_shell_r2.name = 'lower_rod_shell_r2'        # # lower part of the crod with INCO in it - for part with fuel and inco and air only - rod 1
lower_rod_shell_r1.name = 'lower_rod_shell_r1'          # # lower part of the crod with INCO in it - for part with fuel and inco and air only - rod 2

G.name = 'all_graphite_universe'



############################################
### Full 3D core
############################################


##### full core geometry option
geometry = openmc.Geometry([core, upcomer_cell, downcomer_cell, vessel_cell, can_cell,
                            upper_plena_shell_cell, upper_plenum_fuel_cell,
                            lower_plena_shell_cell, lower_plena_fuel_cell,
                            ts_cell, ts_upper_stainless, ts_outer_stainless, ts_hori_cyl_lower_plate, ts_inner_plate, ts_hori_cyl,
                            insulation_cell_upper, insulation_inner_wall, insulation_lower,
                            gas_cell])
geometry.export_to_xml()

############################################
### Partial 3D core (no external)
############################################
# geometry = openmc.Geometry([core, upcomer_cell, downcomer_cell, vessel_cell, can_cell,
#                             upper_plena_shell_cell, upper_plenum_fuel_cell,
#                             lower_plena_shell_cell, lower_plena_fuel_cell])
# geometry.export_to_xml()


############################################
### 2D CORE EZ MODE
############################################
# final_uni = openmc.Universe(cells=[core, upcomer_cell, downcomer_cell, vessel_cell, can_cell,])

# upper_cut = openmc.ZPlane(z0=10, boundary_type = 'reflective')
# lower_cut = openmc.ZPlane(z0=-10, boundary_type = 'reflective')
# core_2d_outer = openmc.ZCylinder(r=outer_vessel_surface.r, boundary_type = 'vacuum')

# final_cell = openmc.Cell(fill=final_uni, region=(-upper_cut & +lower_cut & -core_2d_outer))
# geometry = openmc.Geometry([final_cell])
# geometry.export_to_xml()


############################################
### GEOMETRY PLOTTING
############################################

print("Now plotting some geometry ...")

### PLOTTING YZ
plot = openmc.Plot.from_geometry(geometry)
plot.color_by = 'material'
plot.width = (250,250)
plot.pixels=(500,500)
plot.basis='xz'
plot.origin = (0, L, 0)
plot.colors = {
  f1:'black',
  g1:'pink',
  inconel:'red',
  poison:'green',
  helium:'yellow',
  ss316:'purple',
  inor:'teal',
  tshield:'gray',
  cellgas:'magenta'
}
# plot.to_ipython_image()

### PLOTTING YZ - close up on control rods
plot = openmc.Plot.from_geometry(geometry)
plot.color_by = 'material'
plot.width = (-30,50)
plot.pixels=(5000,5000)
plot.basis='xz'
plot.origin = (0, 0, 0)
plot.colors = {
  f1:'black',
  g1:'pink',
  inconel:'red',
  poison:'green',
  helium:'yellow',
  ss316:'purple',
  inor:'teal',
  tshield:'gray',
  cellgas:'magenta'
}
# plot.to_ipython_image()


### PLOTTING XY
plot = openmc.Plot.from_geometry(geometry)
plot.color_by = 'material'
plot.width = (20, 20)
plot.pixels=[1000,1000]
plot.basis='xy'
# plot.origin = (70, 0, 0)

# plot.colors = {
#   f1:'black',
#   g1:'pink',
#   inconel:'red',
#   poison:'green',
#   helium:'yellow',
#   ss316:'purple',
#   inor:'teal'
# }

#plot.to_ipython_image()


scale252 = [2.000000000E+07	,1.733000000E+07	,1.568000000E+07	,1.455000000E+07	,1.384000000E+07	,
 1.284000000E+07	,1.000000000E+07	,8.187000000E+06	,6.434000000E+06	,4.800000000E+06	,
 4.304000000E+06	,3.679000000E+06	,2.479000000E+06	,2.354000000E+06	,1.850000000E+06	,
 1.500000000E+06	,1.400000000E+06	,1.353000000E+06	,1.317000000E+06	,1.250000000E+06	,
 1.200000000E+06	,1.100000000E+06	,1.010000000E+06	,9.200000000E+05	,9.000000000E+05	,
 8.750000000E+05	,8.611000000E+05	,8.200000000E+05	,7.500000000E+05	,6.790000000E+05	,
 6.700000000E+05	,6.000000000E+05	,5.730000000E+05	,5.500000000E+05	,5.000000000E+05	,
 4.700000000E+05	,4.400000000E+05	,4.200000000E+05	,4.000000000E+05	,3.300000000E+05	,
 2.700000000E+05	,2.000000000E+05	,1.490000000E+05	,1.283000000E+05	,1.000000000E+05	,
 8.500000000E+04	,8.200000000E+04	,7.500000000E+04	,7.300000000E+04	,6.734000000E+04	,
 5.200000000E+04	,5.000000000E+04	,4.500000000E+04	,3.000000000E+04	,2.000000000E+04	,
 1.700000000E+04	,1.300000000E+04	,9.118000000E+03	,8.030000000E+03	,5.700000000E+03	,
 3.900000000E+03	,3.740000000E+03	,3.000000000E+03	,2.500000000E+03	,2.250000000E+03	,
 2.200000000E+03	,1.800000000E+03	,1.550000000E+03	,1.500000000E+03	,1.150000000E+03	,
 9.500000000E+02	,6.830000000E+02	,6.700000000E+02	,5.500000000E+02	,3.050000000E+02	,
 2.850000000E+02	,2.400000000E+02	,2.200000000E+02	,2.095000000E+02	,2.074000000E+02	,
 2.020000000E+02	,1.930000000E+02	,1.915000000E+02	,1.885000000E+02	,1.877000000E+02	,
 1.800000000E+02	,1.700000000E+02	,1.487000000E+02	,1.220000000E+02	,1.190000000E+02	,
 1.175000000E+02	,1.160000000E+02	,1.130000000E+02	,1.080000000E+02	,1.050000000E+02	,
 1.012000000E+02	,9.700000000E+01	,9.000000000E+01	,8.170000000E+01	,8.000000000E+01	,
 7.600000000E+01	,7.200000000E+01	,6.750000000E+01	,6.500000000E+01	,6.300000000E+01	,
 6.100000000E+01	,5.800000000E+01	,5.340000000E+01	,5.060000000E+01	,4.830000000E+01	,
 4.520000000E+01	,4.400000000E+01	,4.240000000E+01	,4.100000000E+01	,3.960000000E+01	,
 3.910000000E+01	,3.800000000E+01	,3.763000000E+01	,3.727000000E+01	,3.713000000E+01	,
 3.700000000E+01	,3.600000000E+01	,3.550000000E+01	,3.500000000E+01	,3.375000000E+01	,
 3.325000000E+01	,3.175000000E+01	,3.125000000E+01	,3.000000000E+01	,2.750000000E+01	,
 2.500000000E+01	,2.250000000E+01	,2.175000000E+01	,2.120000000E+01	,2.050000000E+01	,
 2.000000000E+01	,1.940000000E+01	,1.850000000E+01	,1.700000000E+01	,1.600000000E+01	,
 1.440000000E+01	,1.290000000E+01	,1.190000000E+01	,1.150000000E+01	,1.000000000E+01	,
 9.100000000E+00	,8.100000000E+00	,7.150000000E+00	,7.000000000E+00	,6.875000000E+00	,
 6.750000000E+00	,6.500000000E+00	,6.250000000E+00	,6.000000000E+00	,5.400000000E+00	,
 5.000000000E+00	,4.700000000E+00	,4.000000000E+00	,3.730000000E+00	,3.500000000E+00	,
 3.200000000E+00	,3.100000000E+00	,3.000000000E+00	,2.970000000E+00	,2.870000000E+00	,
 2.770000000E+00	,2.670000000E+00	,2.570000000E+00	,2.470000000E+00	,2.380000000E+00	,
 2.300000000E+00	,2.210000000E+00	,2.120000000E+00	,2.000000000E+00	,1.940000000E+00	,
 1.860000000E+00	,1.770000000E+00	,1.680000000E+00	,1.590000000E+00	,1.500000000E+00	,
 1.450000000E+00	,1.400000000E+00	,1.350000000E+00	,1.300000000E+00	,1.250000000E+00	,
 1.225000000E+00	,1.200000000E+00	,1.175000000E+00	,1.150000000E+00	,1.140000000E+00	,
 1.130000000E+00	,1.120000000E+00	,1.110000000E+00	,1.100000000E+00	,1.090000000E+00	,
 1.080000000E+00	,1.070000000E+00	,1.060000000E+00	,1.050000000E+00	,1.040000000E+00	,
 1.030000000E+00	,1.020000000E+00	,1.010000000E+00	,1.000000000E+00	,9.750000000E-01	,
 9.500000000E-01	,9.250000000E-01	,9.000000000E-01	,8.500000000E-01	,8.000000000E-01	,
 7.500000000E-01	,7.000000000E-01	,6.500000000E-01	,6.250000000E-01	,6.000000000E-01	,
 5.500000000E-01	,5.000000000E-01	,4.500000000E-01	,4.000000000E-01	,3.750000000E-01	,
 3.500000000E-01	,3.250000000E-01	,3.000000000E-01	,2.750000000E-01	,2.500000000E-01	,
 2.250000000E-01	,2.000000000E-01	,1.800000000E-01	,1.400000000E-01	,1.250000000E-01	,
 1.000000000E-01	,9.000000000E-02	,8.000000000E-02	,7.000000000E-02	,6.000000000E-02	,
 5.000000000E-02	,4.200000000E-02	,3.000000000E-02	,2.530000000E-02	,1.000000000E-02	,
 7.500000000E-03	,5.000000000E-03	,4.000000000E-03	,3.000000000E-03	,2.500000000E-03	,
 2.000000000E-03	,1.500000000E-03	,1.200000000E-03	,1.000000000E-03	,7.500000000E-04	,
 5.000000000E-04	,0.000000000E+00]
scale252 = scale252[::-1]



############################################
### CROSS SECTION GENERATION
############################################
# possible group structures from spreadsheet
group_0 =  [0.0000E+00,	6.2500E-01,	2.0000E+07]
group_1 =  [0.0000E+00,	8.0000E-02,	6.2500E-01,	9.1180E+03,	2.0000E+07]
group_2 =  [0.0000E+00,	8.0000E-02,	6.2500E-01,	4.0000E+00,	1.4870E+02,	9.1180E+03,	6.7340E+04,	1.3530E+06,	2.0000E+07]
group_3 =  [0.0000E+00,	4.2000E-02,	8.0000E-02,	4.0000E-01,	6.2500E-01,	4.0000E+00,	1.4870E+02,	9.1180E+03,	6.7340E+04,	5.0000E+05,	1.3530E+06,	3.6790E+06,	2.0000E+07]
group_4 =  [0.0000E+00,	4.2000E-02,	8.0000E-02,	1.4000E-01,	2.5000E-01,	4.0000E-01,	1.3000E+00,	4.0000E+00,	1.4870E+02,	9.1180E+03,	6.7340E+04,	5.0000E+05,	1.3530E+06,	3.6790E+06,	2.0000E+07]
group_5 =  [0.0000E+00,	4.2000E-02,	8.0000E-02,	1.4000E-01,	1.8000E-01,	2.5000E-01,	4.0000E-01,	6.2500E-01,	1.3000E+00,	4.0000E+00,	1.4870E+02,	9.1180E+03,	6.7340E+04,	5.0000E+05,	1.3530E+06,
            3.6790E+06,	2.0000E+07]
group_6 =  [0.0000E+00,	1.2396E-02,	3.5500E-02,	5.6922E-02,	8.1968E-02,	1.1157E-01,	1.4572E-01,	1.8443E-01,	2.2769E-01,	2.5103E-01,	2.9074E-01,	4.5000E-01,	1.4739E+02,	2.0000E+07]
group_7 =  [0.0000E+00,	7.3000E-01,	2.9023E+01,	9.1188E+03,	2.0000E+07]
group_8 =  [0.0000E+00,	1.8554E+00,	2.9023E+01,	9.1188E+03, 2.0000E+07]
group_structures = [group_0, group_1, group_2, group_3, group_4, group_5, group_6, group_7, group_8]
groups = group_structures[user_input_group_struct]
fine_groups = scale252

vols = 0 # dummy
### old mesh
xs_mesh = openmc.RegularMesh()
xs_mesh.dimension = [27,27,xs_gen_axial_intervals]
xs_mesh.lower_left = [-13.5*L, -13.5*L, core_lower_z]
xs_mesh.upper_right = [13.5*L, 13.5*L, core_upper_z]

# ### new mesh
# xs_mesh = openmc.RectilinearMesh()
# # xs_mesh.dimension = [27,27,6]
# xs_mesh.x_grid = np.array([-13.5, -12.5, -11.5, -10.5, -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5])*L
# xs_mesh.y_grid = np.array([-13.5, -12.5, -11.5, -10.5, -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5])*L
# xs_mesh.z_grid = np.array([-85.09000000, -78.10500000, -39.05250000, 0, 39.05250000,  78.10500000, 85.09000000])

num_delayed_groups = 6
legendre = 0
mgxs_types = ['nu-fission', 'nu-scatter', 'reduced absorption', 'scatter probability matrix', # added nu-scatter and reduced absorbtion
  'chi', 'inverse-velocity', 'chi-prompt', 'chi-delayed',  'beta', 'decay-rate', 'kappa-fission']


### list of cells to tally.
cell_domains = [upcomer_cell, downcomer_cell, vessel_cell, can_cell, lower_rod_fuel_r1, lower_rod_fuel_r2, lower_rod_air_r1, lower_rod_air_r2, lower_rod_shell_r1, lower_rod_shell_r2,
                R_only_cell_1, R_only_cell_2, rodded_flibe_1, rodded_flibe_2, upper_plena_shell_cell, upper_plenum_fuel_cell, lower_plena_shell_cell, lower_plena_fuel_cell,
                upperLeft,lowerLeft,upperRight,allButRods]

uni_domains =  [G,
blue_TR, purple_TR, green_TR, yellow_TR, brown_TR,
blue_BR, purple_BR, green_BR, yellow_BR, brown_BR,
blue_BL, purple_BL, green_BL, yellow_BL, brown_BL,
blue_TL, purple_TL, green_TL, yellow_TL, brown_TL] # these tallies for universes - make sure they all have names!!!!
### we want G as our uni domain since this is the one surrounding the core!

tallies_file = openmc.Tallies()
if xs_gen_logical == 1:
  tallies_file, mgxs_mesh_lib, mgxs_mesh_lib_transport, mgxs_cell_lib, mgxs_cell_lib_transport, mgxs_uni_lib, mgxs_uni_lib_transport = make_xs_tallies(xs_mesh, geometry, groups,
                                                                                                                                      mgxs_types, num_delayed_groups,
                                                                                                                                      cell_domains, uni_domains, legendre, fine_groups, tallies_file,
                                                                                                                                      True, True, True)


############################################
### Other tallies
############################################
meshFluxMesh = openmc.RegularMesh()
meshFluxMesh.dimension = [27*1,27*1,50]
meshFluxMesh.lower_left = [-13.5*L, -13.5*L, core_lower_z]
meshFluxMesh.upper_right = [13.5*L, 13.5*L, core_upper_z]
meshFluxMeshFilter = openmc.MeshFilter(meshFluxMesh)
meshFluxEnergyFilter = openmc.EnergyFilter(groups)

flux = openmc.Tally(tally_id=9999, name="fluxMeshTally")
flux.scores = ['flux']
flux.filters.append(meshFluxMeshFilter)
flux.filters.append(meshFluxEnergyFilter)

tallies_file.append(flux)

fissionMesh = openmc.RectilinearMesh(mesh_id=101, name='fission_mesh')
fissionMesh.x_grid = np.array([-13.5*L-rad,     -13.5*L+rad,
                               -13.5*L+1*L-rad, -13.5*L+1*L+rad,
                               -13.5*L+2*L-rad, -13.5*L+2*L+rad,
                               -13.5*L+3*L-rad, -13.5*L+3*L+rad,
                               -13.5*L+4*L-rad, -13.5*L+4*L+rad,
                               -13.5*L+5*L-rad, -13.5*L+5*L+rad,
                               -13.5*L+6*L-rad, -13.5*L+6*L+rad,
                               -13.5*L+7*L-rad, -13.5*L+7*L+rad,
                               -13.5*L+8*L-rad, -13.5*L+8*L+rad,
                               -13.5*L+9*L-rad, -13.5*L+9*L+rad,
                               -13.5*L+10*L-rad, -13.5*L+10*L+rad,
                               -13.5*L+11*L-rad, -13.5*L+11*L+rad,
                               -13.5*L+12*L-rad, -13.5*L+12*L+rad,
                               -13.5*L+13*L-rad, -13.5*L+13*L+rad,
                               -13.5*L+14*L-rad, -13.5*L+14*L+rad,
                               -13.5*L+15*L-rad, -13.5*L+15*L+rad,
                               -13.5*L+16*L-rad, -13.5*L+16*L+rad,
                               -13.5*L+17*L-rad, -13.5*L+17*L+rad,
                               -13.5*L+18*L-rad, -13.5*L+18*L+rad,
                               -13.5*L+19*L-rad, -13.5*L+19*L+rad,
                               -13.5*L+20*L-rad, -13.5*L+20*L+rad,
                               -13.5*L+21*L-rad, -13.5*L+21*L+rad,
                               -13.5*L+22*L-rad, -13.5*L+22*L+rad,
                               -13.5*L+23*L-rad, -13.5*L+23*L+rad,
                               -13.5*L+24*L-rad, -13.5*L+24*L+rad,
                               -13.5*L+25*L-rad, -13.5*L+25*L+rad,
                               -13.5*L+26*L-rad, -13.5*L+26*L+rad,
                               -13.5*L+27*L-rad, -13.5*L+27*L+rad])
fissionMesh.y_grid = fissionMesh.x_grid
fissionMesh.z_grid = np.linspace(core_lower_z, core_upper_z, 51)
fissionMeshFilter = openmc.MeshFilter(fissionMesh)

fission_grid_tally = openmc.Tally(tally_id=100000, name='fission_mesh')
fission_grid_tally.scores = ['fission']
fission_grid_tally.filters.append(fissionMeshFilter)
tallies_file.append(fission_grid_tally)

bypass_cell_filter = openmc.CellFilter([upcomer_cell])
lower_plenum_cell_filter = openmc.CellFilter([lower_plena_fuel_cell])
upper_plenum_cell_filter = openmc.CellFilter([upper_plenum_fuel_cell])
downcomer_cell_filter = openmc.CellFilter([downcomer_cell])

fissions_downcomer = openmc.Tally(tally_id=100001, name='fissions_downcomer')
fissions_UP = openmc.Tally(tally_id=100002, name='fissions_UP')
fissions_LP = openmc.Tally(tally_id=100003, name='fissions_LP')
fissions_bypass = openmc.Tally(tally_id=100004, name='fissions_bypass')

fissions_downcomer.scores = ['fission']
fissions_UP.scores = ['fission']
fissions_LP.scores = ['fission']
fissions_bypass.scores = ['fission']

fissions_downcomer.filters.append(downcomer_cell_filter)
fissions_UP.filters.append(upper_plenum_cell_filter)
fissions_LP.filters.append(lower_plenum_cell_filter)
fissions_bypass.filters.append(bypass_cell_filter)

tallies_file.append(fissions_downcomer)
tallies_file.append(fissions_LP)
tallies_file.append(fissions_UP)
tallies_file.append(fissions_bypass)

tallies_file.export_to_xml()
############################################
### VOLUME CALCULATION
############################################
lower_left_lowplena =     (-32*2.54, -32*2.54, -67*2.54/2-50)
upper_right_lowplena =    (32*2.54, 32*2.54, -67*2.54/2+2)
lower_left_upperplena =   (-32*2.54, -32*2.54, 67*2.54/2-2)
upper_right_upperplena =  (32*2.54, 32*2.54, 67*2.54/2+50)

num_vol_samples = 5000000000 # 5 billion neutrons - approximate runtime == 22.81*50 = 1000 seconds = 19 minutes on 5 cores

vol_calc_lower_plena = openmc.VolumeCalculation([ lower_plena_fuel_cell], num_vol_samples, lower_left_lowplena, upper_right_lowplena)
vol_calc_upper_plena = openmc.VolumeCalculation([ upper_plenum_fuel_cell], num_vol_samples, lower_left_upperplena, upper_right_upperplena)


############################################
### SETTINGS
############################################
source = openmc.Source()
source.space = openmc.stats.Box((-70, -70, -170/2), (70, 70, 170/2))

settings = openmc.Settings()
settings.volume_calculations = [vol_calc_lower_plena, vol_calc_upper_plena]
# settings.run_mode = 'volume' # turn on for vol calcs
settings.source = source
settings.batches = nact+nsk
settings.inactive = nsk
settings.particles = npg
settings.output = {'tallies': False}
settings.temperature = {"method": "interpolation"}
settings.export_to_xml()
# !cat settings.xml


############################################
### RUN OPENMC
############################################
# openmc.run(threads=int(nthreads))


############################################
###
############################################


############################################
###
############################################
