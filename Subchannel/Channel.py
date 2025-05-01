import numpy as np
from .FluidRelation import FluidRelation
from Meshing.Meshing import *
from Fields.Fields import *
from Kernels.Kernels import *
from Solvers.Solvers import *

class Channel:
  def __init__(self, gravity: float,
               Dh: float,
               area: float,
               temp_tolerance: float,
               max_temp_iterations: int,
               nZones: int,
               L0: float,
               L1: float,
               fluid: FluidRelation,
               pressure_bc: float,
               T_bc: float,
               mdot_bc: float,
               fric: str,
               heat_source: list | float):
    # Initial data that remains constant for this channel during the simulation.
    self.gravity = gravity
    self.Dh = Dh
    self.area = area
    self.temp_tolerance = temp_tolerance
    self.max_temp_iterations = max_temp_iterations
    self.nZones = nZones
    self.L0 = L0
    self.L1 = L1
    self.fluid = fluid
    self.fric = fric
    self.set_heat_source(nZones=nZones, heat_source=heat_source)

    # Make a mesh
    coords = np.linspace(L0,L1,self.nZones+1)
    areas = [area]*(self.nZones+1)
    self.mesh = Mesh_1D(nodeCoords=coords, faceAreas=
                        areas)
    # Compute total channel volume
    self.ch_volume = 0.0
    for cid in self.mesh.cidList:
      self.ch_volume += self.mesh.cells[cid].vol

    # Set boundary conditions
    self.set_bcs(pressure_bc=pressure_bc, T_bc=T_bc, mdot_bc=mdot_bc, tracer_name_value_pairs={}, tracer_bool=False, th_bool=True)

    # Make some variables
    self.pressure = ScalarField(name='P', initial_value=self.pressure_bc, mesh=self.mesh)
    self.mdot = ScalarField(name='mdot', initial_value=self.mdot_bc, mesh=self.mesh)
    self.h = ScalarField(name='enthalpy', initial_value=self.T_bc, mesh=self.mesh)
    self.rho = ScalarField(name='density', initial_value=self.rho_bc, mesh=self.mesh)
    self.temp = ScalarField(name='temp', initial_value=self.T_bc, mesh=self.mesh)

    # Make a dictionary for those variables.
    self.pressure_dict = {}
    self.mdot_dict = {}
    self.h_dict = {}
    self.rho_dict = {}
    self.temp_dict = {}

    # Make some matrices and vectors
    self.A_mass = np.zeros([self.mesh.nz, self.mesh.nz])
    self.A_energy = np.zeros([self.mesh.nz, self.mesh.nz])
    self.A_momentum = np.zeros([self.mesh.nz, self.mesh.nz])
    self.b_mass = np.zeros([self.mesh.nz])
    self.b_energy = np.zeros([self.mesh.nz])
    self.b_momentum = np.zeros([self.mesh.nz])

    self.A_pressure_mass = np.zeros([self.mesh.nz*2, self.mesh.nz*2])
    self.b_pressure_mass = np.zeros([self.mesh.nz*2])

    # Channel conditions dictionary for passing data between channels
    self.channel_conditions = {}
    self.channel_conditions['dP'] = None
    self.channel_conditions['P_out'] = None
    self.channel_conditions['P_in'] = None
    self.channel_conditions['h_out'] = None
    self.channel_conditions['h_in'] = None
    self.channel_conditions['mdot_out'] = None
    self.channel_conditions['mdot_in'] = None
    self.channel_conditions['T_in'] = None
    self.channel_conditions['T_out'] = None


    # Channel face field handling -> FaceField for velocity essentially.
    self.velocity_faces = FaceField(name='vel', initial_value=0.0 , mesh=self.mesh) # nZones + 1

    # Channel Tracers -> keys are tracer names
    self.tracers = {} # Field variable for tracers
    self.tracer_kernels = {} # actual kernel objects for tracers
    self.tracer_bcs = {} # actual BC objects for tracers
    self.tracer_bc_values = {} # name and value pairs for tracers [name] -> inlet bc value for this channel
    self.channel_conditions['tracers_in'] = {} # {} -> name and value pairs for inlet values and outlet values
    self.channel_conditions['tracers_out'] = {} #
    self.tracer_dict = {} # time dependent tracer dict -> tracer_dict[tracer_name][_t] -> vector of values for this tracer


  # SETUP A TRACER FOR THIS CHANNEL OBJECT
  def add_tracer_to_channel(self, name: str,
                 initial_value: np.ndarray | float,
                 scheme: str,
                 decay_const: float,
                 boundary: str,
                 phi: float,
                 rho: float,
                 source: float | np.ndarray | ScalarField,
                 beta: float):
    """
    Adds a tracer to the channel - advection + decay + source
    """

    # Raise Exception (only 1 way advection supported currently since I am lazy and it really doesnt matter for 1D)
    if boundary != 'lower':
      raise Exception("Only lower boundary conditions supported (inlet at bottom and zero gradient outlet at top!)")

    # Add a tracer object
    self.tracers[name] = ScalarField(name=name, initial_value=initial_value, mesh=self.mesh)

    # Add kernels to the tracer:
    self.tracer_kernels[name] = [AdvectionKernel(field=self.tracers[name], mesh=self.mesh, w=self.velocity_faces, scheme=scheme, rho=rho),
                                 ImplicitReactionKernel(field=self.tracers[name], mesh=self.mesh, lam=decay_const),
                                 ExplicitSourceKernel(field=self.tracers[name], mesh=self.mesh, Q=source, beta=beta)]

    # Add BC's to the tracer:
    self.tracer_bcs[name] = [AdvectedInletFluxBC(field=self.tracers[name], mesh=self.mesh, boundary=boundary, phi=phi, w=self.velocity_faces, rho=rho)]

  # SOLUTION FOR CHANNEL TRACER
  def solve_tracer(self, name: str, _dt: float):
    """
    Solves tracer equations.
    """
    # Setup basic solver
    solver = BasicSolver(kernels=self.tracer_kernels[name], bcs=self.tracer_bcs[name], field=self.tracers[name])

    # Solve
    solver.solve()

    # Update channel conditions dictionary at the inlet and outlet.
    self.channel_conditions['tracers_in'][name] = self.tracer_bcs[name][0].phi # inlet value from BC
    self.channel_conditions['tracers_out'][name] = self.tracers[name].T[-1] # outlet value using zero gradient at the outlet

  def solve_all_tracers(self, _dt: float):
    for tracer_name in self.tracers.keys():
      self.solve_tracer(name=tracer_name, _dt=_dt)

  # CONVERT MASS FLOW RATES TO VELOCITY FACEFIELD
  def mdot_to_velocities(self):
    """
    Def. that passes takes current values of mdot and updates velocity faces field
    """
    # Assume inlet bc is known - face 0
    # Assume zero gradient at outlet - face -1
    self.velocity_faces.T[0] = self.mdot_bc / self.rho_bc / self.area
    for cid in self.mesh.cidList:
      if cid == 0:
        continue
      else:
        gC = self.mesh.cells[cid-1].geoUpper
        mdot_interp = self.mdot.T[cid-1] * gC + (1.0-gC)*self.mdot.T[cid]
        rho_interp = self.rho.T[cid-1] * gC + (1.0-gC)*self.rho.T[cid]
        face_area = self.area
        vel_f = mdot_interp/rho_interp/face_area
        self.velocity_faces.T[cid] = vel_f
    mdot_last = self.mdot.T[-1]
    rho_last = self.rho.T[-1]
    area_last = self.area
    self.velocity_faces.T[-1] = mdot_last/rho_last/area_last

  # UPDATE OLD VARIABLE TO NEW VALUE
  def update_old_to_most_recent(self):
    """
    Updates old -> new
    """
    # THERMAL HYDRAULIC VARIABLES
    self.pressure.T_old = copy.deepcopy(self.pressure.T)
    self.mdot.T_old = copy.deepcopy(self.mdot.T)
    self.h.T_old = copy.deepcopy(self.h.T)
    self.rho.T_old = copy.deepcopy(self.rho.T)
    self.temp.T_old = copy.deepcopy(self.temp.T)

    # TRACER VARIABLES
    for name in self.tracers.keys():
      self.tracers[name].T_old = copy.deepcopy(self.tracers[name].T)

  # SAVES DEEPCOPY OF DATA FOR STORAGE/USE LATER
  def save_data(self, _t: float):
    """
    Stores data for the next timestep - called externally.
    """

    # Saves data
    self.pressure_dict[_t] = copy.deepcopy(self.pressure.T)
    self.mdot_dict[_t] = copy.deepcopy(self.mdot.T)
    self.h_dict[_t] = copy.deepcopy(self.h.T)
    self.rho_dict[_t] = copy.deepcopy(self.rho.T)
    self.temp_dict[_t] = copy.deepcopy(self.temp.T)

    # Saves data for tracers
    for name in self.tracers.keys():
      self.tracer_dict[name][_t] = copy.deepcopy(self.tracers[name].T)

  # SETS HEAT SOURCE
  def set_heat_source(self, heat_source: float | list, nZones: int):
    # Quick function for setting heat source.
    if isinstance(heat_source, float):
      self.heat_source = [heat_source]*nZones # float
    else:
      self.heat_source = heat_source # if heat source is a list

  # FRICTION FACTOR CORRELATIONS
  def get_friction_factor(self, Reynolds: float):
    # Friction  factors for a channel from the thesis Luzzi et al., 2010
    if self.fric == 'type1':
      if Reynolds > 3000:
        return 0.3164 / Reynolds**0.25 * (1.0 + Reynolds/4.31e5)**(1.0/8.0)
      elif Reynolds <= 3000:
        return 64.0 / Reynolds
  # CASE OF NO FRICTION FACTOR CORRELATIONS
    if self.fric == 'none':
      return 0.0
    # Unknown!
    else:
      raise Exception("Friction factor type unknown!")

  # UPDATE OR SET BOUNDARY CONDITIONS
  def set_bcs(self, pressure_bc: float, T_bc: float, mdot_bc: float, tracer_name_value_pairs: dict, tracer_bool: bool, th_bool: bool):
    """
    Sets boundary conditions for the pressure, temperature, and mdot.
    """
    # TH BOUNDARY CONDITIONS
    if th_bool:
      self.pressure_bc = pressure_bc
      self.T_bc = T_bc
      self.mdot_bc = mdot_bc
      self.h_bc = self.fluid.props_from_P_T(P=self.pressure_bc, T=self.T_bc, prop='h')
      self.rho_bc = self.fluid.props_from_P_H(P=self.pressure_bc, enthalpy=self.h_bc, prop='rho')

    # TRACER BOUNDARY CONDITIONS - incoming values for C_i for the tracer.
    if tracer_bool:
      self.tracer_bc_values = tracer_name_value_pairs # key -> tracer name , value -> bc value for that tracer
      for tracer_name in tracer_name_value_pairs.keys():
        self.tracer_bcs[tracer_name][0].phi = tracer_name_value_pairs[tracer_name] # use index 0 since bcs is a list of bc's

  # MASS EQUATION
  def solve_mass_equation(self, _dt: float):
    """
    MASS EQUATION AND SETUP
    """
    mesh = self.mesh
    self.A_mass *= 0.0
    self.b_mass *= 0.0
    for cid in mesh.cidList:
      area = ( mesh.cells[cid].upperArea + mesh.cells[cid].lowerArea ) / 2.0
      dz = mesh.cells[cid].dz
      A_X_dt = area*dz / _dt

      # Solve implicitely for mass flow rates based on densities
      ### A_X_dt * (rho - rho_old) + (m - m_i-1) = 0.0
      self.b_mass[cid] =  -1.0 * A_X_dt * (self.rho.T[cid] - self.rho.T_old[cid])
      self.A_mass[cid,cid] = 1.0
      if cid == 0: # Boundary
        self.b_mass[cid] += self.mdot_bc
      else:
        self.A_mass[cid,cid-1] = -1.0

    # SOLVE
    self.mdot.T = np.linalg.solve(self.A_mass, self.b_mass)

  # ENERGY EQUATION
  def solve_energy_equation(self, _dt: float):
    """
    ENERGY EQUATION AND SETUP
    """
    mesh = self.mesh
    self.A_energy *= 0.0
    self.b_energy *= 0.0

    for cid in mesh.cidList:
      # Geometry stuff
      area = ( mesh.cells[cid].upperArea + mesh.cells[cid].lowerArea ) / 2.0
      dz = mesh.cells[cid].dz
      # main diagonal
      self.A_energy[cid,cid] = area/_dt * self.rho.T[cid] + 1.0/dz * self.mdot.T[cid]

      # Source term (heat source in W/m3 --> times Volume ---> divided by dz ---> = linear heat rate)
      self.b_energy[cid] += self.heat_source[cid] * area*dz / dz

      # Time term for source
      self.b_energy[cid] += area/_dt * self.h.T_old[cid] * self.rho.T_old[cid]

      # off diagonal coeff for enthalpy
      if cid == 0:
        self.b_energy[cid] += 1.0/dz * self.mdot_bc * self.h_bc
      else:
        self.A_energy[cid,cid-1] = -1.0/dz * self.mdot.T[cid-1]

    # SOLVE
    self.h.T = np.linalg.solve(self.A_energy, self.b_energy)

  # PRESSURE EQUATION
  def solve_pressure_equation(self, _dt: float):
    """
    MOMENTUM EQUATION AND SETUP
    """
    mesh = self.mesh
    self.A_momentum *= 0.0
    self.b_momentum *= 0.0
    for cid in mesh.cidList:
      # geometry stuff
      area = ( mesh.cells[cid].upperArea + mesh.cells[cid].lowerArea ) / 2.0
      dz = mesh.cells[cid].dz

      # main diagonal
      self.A_momentum[cid, cid] = area

      # off diagonal
      if cid == 0:
        self.b_momentum[cid] = self.pressure_bc*area
      else:
        self.A_momentum[cid,cid-1] = -1.0*area

      # compute friction factor
      Reynolds = self.mdot.T[cid] * self.Dh * self.rho.T[cid] / self.fluid.get_mu()
      _fric = self.get_friction_factor(Reynolds=Reynolds)

      # additional b terms
      self.b_momentum[cid] += (
            -1.0 * self.gravity*area*dz*self.rho.T[cid]) - (
            0.5 * dz * _fric/self.Dh/self.rho.T[cid] * np.abs(self.mdot.T[cid])  * self.mdot.T[cid] / area ) - (
            dz / _dt * (self.mdot.T[cid] - self.mdot.T_old[cid]) )
    self.pressure.T = np.linalg.solve(self.A_momentum, self.b_momentum)

  # EQUATION OF STATE FOR TEMPERATURE AND DENSITY
  def solve_EOS(self):
    # New temp and density arrays
    temp_updated = np.array([])
    rho_updated = np.array([])
    for cid in self.mesh.cidList:
      t_new = self.fluid.props_from_P_H(P=self.pressure.T[cid], enthalpy=self.h.T[cid], prop='T')
      temp_updated = np.append(temp_updated, t_new)
      rho_updated = np.append(rho_updated, self.fluid.props_from_P_H(P=self.pressure.T[cid], enthalpy=self.h.T[cid], prop='rho') )

    max_diff = np.max(np.abs(self.temp.T - temp_updated))
    self.temp.T = temp_updated
    self.rho.T = rho_updated
    return max_diff

  # UPDATES DICT FOR INLET AND OUTLET CONDITIONS.
  def update_channel_conditions_TH(self):
    """
    Called at the end of a TH loop to update TH incoming and outgoing fluxes/values.
    """
    self.channel_conditions['P_out'] = self.pressure.T[-1]
    self.channel_conditions['P_in'] = self.pressure_bc
    self.channel_conditions['dP'] = self.pressure.T[-1] - self.pressure_bc
    self.channel_conditions['h_out'] = self.h.T[-1]
    self.channel_conditions['h_in'] = self.h_bc
    self.channel_conditions['mdot_out'] = self.mdot.T[-1]
    self.channel_conditions['mdot_in'] = self.mdot_bc
    self.channel_conditions['T_out'] = self.temp.T[-1]
    self.channel_conditions['T_in'] = self.T_bc

  # THERMAL HYDRAULICS SOLVING LOOP
  def solve_channel_TH(self, _dt: float):
    # Start loop
    iteration_num = 0
    max_diff = 1e10
    while max_diff > self.temp_tolerance:

      # Handle iteration limit
      iteration_num += 1
      if iteration_num > self.max_temp_iterations:
        raise Exception("Iteration number exceeded for this channel!")

      # Solve equations
      self.solve_mass_equation(_dt=_dt) # solves and updates mass flow rates
      self.solve_energy_equation(_dt=_dt) # solves and updates enthalpy
      self.solve_pressure_equation(_dt=_dt) # solves and updates pressures
      max_diff = self.solve_EOS() # solves and updates temperature and density - returns maximum temp diff.

    # Print?
    print("Channel solved after", iteration_num, "iterations!")

    # After we solve we update our channel conditions for later reference:
    self.update_channel_conditions_TH()

    # Now update velocity face-fields from mdot and density solution
    self.mdot_to_velocities()


class ChannelArray:
  """
    Class of coupled channels in a channel array.
  """
  def __init__(self, channels: np.ndarray, coupling_method: str, flow_ratios, fluid: FluidRelation):
    # channels is a np array of subchannels
    self.channels = channels

    # coupling method either prescribed mass flow ratios for each channel OR pressure based
    self.coupling_method = coupling_method
    if (self.coupling_method != 'pressure_method') & (self.coupling_method != 'ratio_method'):
      raise Exception("Unknown method for coupling channel array!")

    # Boundary condition at the inlet of the channel array
    self.mdot_bc = None # TOTAL mass flow rate into all channels
    self.pressure_bc = None
    self.T_bc = None
    self.h_bc = None
    self.rho_bc = None

    # Inlet BC dictionary for tracers:
    self.tracer_bcs = {}

  def set_bcs(self, pressure_bc: float, T_bc: float, mdot_bc: float, tracer_name_value_pairs: dict, tracer_bool: bool, th_bool: bool):
    """
    Called by ChannelInterface object to update/bcs for a channel.
    """
    # HANDLING THERMAL HYDRAULIC BOUNDARY CONDITIONS
    if th_bool:
      self.mdot_bc = mdot_bc # TOTAL mass flow rate across all coupled channels
      self.pressure_bc = pressure_bc # pressure value
      self.T_bc = T_bc # advected temperature value
      self.h_bc = self.fluid.props_from_P_T(P=self.pressure_bc, T=self.T_bc, prop='h')
      self.rho_bc = self.fluid.props_from_P_H(P=self.pressure_bc, enthalpy=self.h_bc, prop='rho')

    # HANDLING TRACER BOUNDARY CONDITIONS
    if tracer_bool:
      self.tracer_bcs = copy.deepcopy(tracer_name_value_pairs)
      for key in self.tracer_bcs:
        for ch in self.channels:
          ch.tracer_bcs[key].phi = self.tracer_bcs[key]

  def get_outlet_conditions(self):
    """
    Gets outlet conditions for the array of channels
    """
    advected_enthalpy = 0.0
    pressure_sum = 0.0
    mdot_sum = 0.0
    Nchannels = len(self.channels)

    for ch in self.channels:
      pressure_sum += ch.pressure.T[-1]
      advected_enthalpy += ch.h.T[-1] * ch.mdot.T[-1]
      mdot_sum += ch.mdot.T[-1]

    # outlet pressure - averaged
    outlet_P = pressure_sum/Nchannels

    # advected enthalpy to specific enthalpy
    specific_enthalpy = advected_enthalpy / mdot_sum

    # then get outlet temperature
    outlet_T = self.fluid.props_from_P_H(P=outlet_P, enthalpy=specific_enthalpy, prop='T')

    # Now do advected tracers:
    tracer_outgoing_advected = {}
    tracer_weighted_value = {} # C_out = int(C*A*U)_out / int(A*U)_out
    for tracer_key in self.channels[0].tracers.keys():
      advected_tracer = 0.0
      area_total = 0.0
      vol_flux = 0.0
      for ch in self.channels:
        # basic stuffs for this channel
        vel = ch.velocity_faces.T[-1]
        A = ch.area

        # compute toitals
        vol_flux += vel * ch.area
        area_total += ch.area
        advected_tracer_quantity += ch.tracers[tracer_key].T[-1] * vel * A

      tracer_outgoing_advected[tracer_key] = advected_tracer_quantity
      tracer_weighted_value[tracer_key] = advected_tracer_quantity / vol_flux


    return specific_enthalpy, outlet_P, outlet_T, mdot_sum, tracer_weighted_value

  def solve_channel_TH(self, _dt: float):
    """
    Solves thermal hydraulics for each channel in the array.
    TODO need to include mass flow coupling/iterations.
    TODO need to properly handle bc's in solving channel TH --- mostly done in set_bcs where
         TODO we have to properly setup boundary conditions for each individual channel before we solve it here.
    """
    for ch in self.channels:
      ch.solve_channel_TH(_dt=_dt)

  def solve_tracer(self, name: str, _dt: float):
    """
    Solves Tracers for each channel in the array
    """
    for ch in self.channels:
      ch.solve_tracer(name=name, _dt=_dt)

  def add_tracer_to_channel(self, name: str,
                 initial_value: np.ndarray | float,
                 scheme: str,
                 decay_const: float,
                 boundary: str,
                 phi: float,
                 rho: float,
                 source: float | np.ndarray | ScalarField,
                 beta: float):
    """
    Sets up tracers for all channels in this channel array.
    Do not do this for transient cases since initial conditions are likely not the same for every single channel. Instead
    use a pre-existing ChannelArray object
    """
    for ch in self.channels:
      ch.add_tracer_to_channel(name=name, initial_value=initial_value,
                               scheme=scheme,decay_const=decay_const,
                               boundary=boundary,phi=phi,rho=rho,source=source, beta=beta)

class ChannelInterface:
  """
  Class that serves as an interface between two channels.

  Going to try to generalize things but in general lets just say that flow goes from ch1 to ch2

  note: Index 0 is always the inlet and index -1 is always the outlet.

  Basically just passes data from the output of a channel type to the input of another.

  Ch1 is the outputting channel.
  Ch2 is the receiving channel.
  """
  def __init__(self, ch1, ch2):
    # Channels
    self.ch1 = ch1
    self.ch2 = ch2

  def update_interface_conditions(self, tracer_bool: bool, th_bool: bool):
    """
    Sets outgoing values of main channel (ch1) to incoming values of other channel (ch2)
    """
    # Channel to channel scenario
    if isinstance(self.ch1, Channel) & isinstance(self.ch2, Channel):
      cond1 = self.ch1.channel_conditions
      self.ch2.set_bcs(pressure_bc=cond1['P_out'], T_bc=cond1['T_out'], mdot_bc=cond1['mdot_out'], tracer_name_value_pairs=cond1['tracers_out'],
                       tracer_bool=tracer_bool, th_bool=th_bool)

    # Channel to channel array scenario
    elif isinstance(self.ch1, Channel) & isinstance(self.ch2, ChannelArray):
      # This means that Channel1 passes data to ChannelArray2
      # E.g. the case of lower plenum to core subchannels.
      cond1 = self.ch1.channel_conditions
      self.ch2.set_bcs(pressure_bc=cond1['P_out'], T_bc=cond1['T_out'], mdot_bc=cond1['mdot_out'], tracer_name_value_pairs=cond1['tracers_out'],
                       tracer_bool=tracer_bool, th_bool=th_bool)

    # Channel array to channel scenario
    elif isinstance(self.ch1, ChannelArray) & isinstance(self.ch2, Channel):
      _, outlet_P, outlet_T, mdot_sum, tracer_weighted_outlet_values = self.ch1.get_outlet_conditions()
      self.ch2.set_bcs(pressure_bc=outlet_P, T_bc=outlet_T, mdot_bc=mdot_sum, tracer_name_value_pairs=tracer_weighted_outlet_values,
                       tracer_bool=tracer_bool, th_bool=th_bool)
    else:
      raise Exception("Unknown interface type")

