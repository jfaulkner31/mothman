"""
Imports
"""
import numpy as np
from Kernels.TimeDelayTracker import *

"""
HX object
"""

class HeatExchangerObject:
  def __init__(self, area: float,
               htc: float,
               hot_dict: dict,
               cold_dict: dict,
               temp_rel_tol: float,
               tracer_hot_delay: float = -1.0,
               tracer_cold_delay: float = -1.0):
    """
    Labelling
    ==========
    cd : hot side
    cl : cold side
    """
    self.Tcd_in = None
    self.Tcd_out = None
    self.Tcl_in = None
    self.Tcl_out = None

    self.T_cd_in_old = None
    self.T_cl_in_old = None
    self.T_cd_out_old = None
    self.T_cl_out_old = None

    # Need to do how to handle inputs.
    # especially the input from the main loop.

    # Need enthalpy converter as a Special type of linker
    # for the input

    # Then need time delayers and linkers for others...

    """
    Dicts for hot and cold sides:
    mdot : mass flow rate
    cp : specific heat
    M : mass
    pressure : pressure
    """
    self.hot: dict = hot_dict
    self.cold: dict = cold_dict

    self.A: float = area
    self.htc: float = htc

    self.rel_tol = temp_rel_tol

    """Tracer related"""
    self.hot_tracer_delay: float = tracer_hot_delay
    self.cold_tracer_delay: float = tracer_cold_delay
    self.hot_delayers: list[ExponentialDelayTracker] = []
    self.cold_delayers: list[ExponentialDelayTracker] = []

  def print_solution(self):
    print(f"T_hot_out={self.Tcd_out} | T_cold_out={self.Tcl_out}")
    print(f"T_hot_in={self.Tcd_in} | T_cold_in={self.Tcl_in}")
    print(f"dT_hot={self.Tcd_out-self.Tcd_in} | dT_cold={self.Tcl_out-self.Tcl_in}")
    print(f"LMTD={self._logmean()}")

  def get_hot_outlet_conds(self) -> tuple[float, float, float]:
    """
    Gets hot outlet conditions

    Returns
    =======
    pressure : float
      pressure
    mdot : float
      mass flow rate
    T : float
      temperature of the outlet
    """
    return self.hot['pressure'], self.hot['mdot'], self.Tcd_out

  def get_cold_outlet_conds(self) -> tuple[float, float, float]:
    """
    Gets cold outlet conditions

    Returns
    =======
    pressure : float
      pressure
    mdot : float
      mass flow rate
    T : float
      temperature of the outlet
    """
    return self.cold['pressure'], self.cold['mdot'], self.Tcl_out

  def initialize_outlets(self, transferred_heat: float):
    # Hot side: Q/mdot*cp = dT = Tin - Tout
    self.Tcd_out = self.Tcd_in - transferred_heat/self.hot['mdot']/self.hot['cp']
    # Cold side
    self.Tcl_out = self.Tcl_in + transferred_heat/self.cold['mdot']/self.cold['cp']

    # Set the old guys
    self.Tcl_out_old = self.Tcl_out
    self.T_cd_out_old = self.Tcd_out

  def set_hot_bc(self, T_in: float, mdot: float, pressure: float):
    if self.T_cd_in_old == None:
      self.T_cd_in_old = T_in
    self.Tcd_in = T_in
    self.hot['mdot'] = mdot
    self.hot['pressure'] = pressure

  def set_cold_bc(self, T_in: float, mdot: float, pressure: float):
    if self.T_cl_in_old == None:
      self.T_cl_in_old = T_in
    self.Tcl_in = T_in
    self.cold['mdot'] = mdot
    self.cold['pressure'] = pressure

  def update_old_to_most_recent(self):
    """
    Called externally, updates old to most recent.
    Usually called when soln. is done.
    """
    self.T_cd_in_old = self.Tcd_in
    self.T_cl_in_old = self.Tcl_in
    self.T_cd_out_old = self.Tcd_out
    self.T_cl_out_old = self.Tcl_out

  def solve(self, _dt: float, printSolve: bool = False):
    """
    Solves heat exchanger

    Params
    ======
    _dt : float
      time step in seconds
    """
    lmtd_old = self._logmean()
    diff = 1e6
    it = 0
    while diff > self.rel_tol:
      it += 1
      T_cd_out = self._hot_diffeq(_dt=_dt, info=self.hot)
      T_cl_out = self._cold_diffeq(_dt=_dt, info=self.cold)
      self.Tcd_out = T_cd_out
      self.Tcl_out = T_cl_out
      lmtd_new = self._logmean()
      diff = np.abs(1.0 - lmtd_old/lmtd_new)

      if printSolve:
        print(f"Iteration {it} | diff = {diff} | lmtd_old = {lmtd_old} | lmtd_new = {lmtd_new}")

      lmtd_old = lmtd_new

  def _logmean(self) -> float:
    """
    Returns logmean temperature difference based on
    current temperature values.

    Params
    ======
    None
    """
    top = self.Tcd_in - self.Tcl_in - (self.Tcd_out - self.Tcl_out)
    bottom = np.log(self.Tcd_in - self.Tcl_in) - np.log(self.Tcd_out - self.Tcl_out)
    return top / bottom

  def _heat_source_term(self, info: dict):
    M =  info['mass']
    cp = info['cp']
    return 2.0 / M / cp * self._get_H() * self.A * self.htc * self._logmean()

  def get_heat_transferred(self):
    return self.A * self._get_H() * self.htc * self._logmean()

  def _get_H(self):
    """
    Returns heat exchanger factor based on design
    and current solution.
    """
    return 1.0

  def _hot_diffeq(self, _dt: float, info: dict) -> float:
    """
    Differential equation for the hotside.
    """
    M =  info['mass']
    mdot = info['mdot']
    cp = info['cp']

    F = 1.0/_dt
    C1 = 2.0 * mdot / M
    C2 = self._heat_source_term(info=info)

    Tdiff = (self.Tcd_in - self.T_cd_in_old) * F

    T_cd_out = (-C1 * self.Tcd_in + C2 - Tdiff) / (F - C1)
    return T_cd_out

  def _cold_diffeq(self, _dt: float, info: dict) -> float:
    """
    Differential equation for the cold side
    """
    M =  info['mass']
    mdot = info['mdot']
    cp = info['cp']

    F = 1.0/_dt
    C1 = 2.0 * mdot / M
    C2 = self._heat_source_term(info=info)

    Tdiff = (self.Tcl_in - self.T_cl_in_old) * F

    T_cl_out = (-C1 * self.Tcl_in - C2 - Tdiff) / (F - C1)
    return T_cl_out

  """
  Tracer related methods for HeatExchanger
  """

  """
  Setting bc values for tracers
  """
  def set_tracer_hot_bc(self, current_time: float, tracer_name_value_pairs: dict):
    """
    Sets inlet value at the hot side for the tracers

    Params
    ======
    current_time : float
      the current simulation time in seconds
    tracer_name_value_pairs : dict[str, float]
      dict with names of tracers with values being the value to assign
    """
    for name in tracer_name_value_pairs:
      self.hot_delayers[name].add_pair(current_time, tracer_name_value_pairs[name])

  def set_tracer_cold_bc(self, current_time: float, tracer_name_value_pairs: dict):
    """
    Sets inlet value at the cold side for the tracers

    Params
    ======
    current_time : float
      the current simulation time in seconds
    tracer_name_value_pairs : dict[str, float]
      dict with names of tracers with values being the value to assign
    """
    for name in tracer_name_value_pairs:
      self.cold_delayers[name].add_pair(current_time, tracer_name_value_pairs[name])


  """
  Initial setup of hotside and coldside tracers
  """
  def setup_tracers_hotside(self, reference_channel):
    """
    Sets up tracers for heat exchanger hot side

    Params
    ======
    reference_channel : Channel
      reference channel to get tracer information from
    """
    if self.hot_tracer_delay < 0.0:
      raise Exception("Hot tracer delay not set")
    delayers: dict[str, ExponentialDelayTracker] = {}
    for name in reference_channel.tracer_kernels.keys():
      decay = reference_channel.tracer_kernels[name][1].lam
      delayers[name] = ExponentialDelayTracker(decay=decay,
                                               delay=self.hot_tracer_delay)
    self.hot_delayers = delayers

  def setup_tracers_coldside(self, reference_channel):
    """
    Sets up tracers for heat exchanger cold side

    Params
    ======
    reference_channel : Channel
      reference channel to get tracer information from
    """
    if self.cold_tracer_delay < 0.0:
      raise Exception("Cold tracer delay not set")
    delayers: dict[str, ExponentialDelayTracker] = {}
    for name in reference_channel.tracer_kernels.keys():
      decay = reference_channel.tracer_kernels[name][1].lam
      delayers[name] = ExponentialDelayTracker(decay=decay,
                                               delay=self.cold_tracer_delay)
    self.cold_delayers = delayers

  def _get_hot_tracer_outlet(self, name: str, time: float) -> float:
    """
    Gets hot tracer outlet value

    Params
    ======
    name : str
      tracer name
    time : float
      current time in seconds

    Returns
    =======
    value : float
      tracer value at outlet
    """
    return self.hot_delayers[name].get_delayed_value(time)

  def get_hot_tracer_all(self, time: float) -> dict[str, float]:
    """
    Gets all hot tracer outlet values

    Params
    ======
    time : float
      current time in seconds

    Returns
    =======
    values : dict[str, float]
      tracer name to tracer value at outlet
    """
    values: dict[str, float] = {}
    for name in self.hot_delayers.keys():
      values[name] = self._get_hot_tracer_outlet(name, time)
    return values

  def get_cold_tracer_all(self, time: float) -> dict[str, float]:
    """
    Gets all cold tracer outlet values

    Params
    ======
    time : float
      current time in seconds

    Returns
    =======
    values : dict[str, float]
      tracer name to tracer value at outlet
    """
    values: dict[str, float] = {}
    for name in self.cold_delayers.keys():
      values[name] = self._get_cold_tracer_outlet(name, time)
    return values

  def _get_cold_tracer_outlet(self, name: str, time: float) -> float:
    """
    Gets cold tracer outlet value
    Params
    ======
    name : str
      tracer name
    time : float
      current time in seconds

    Returns
    =======
    value : float
      tracer value at outlet
    """
    return self.cold_delayers[name].get_delayed_value(time)

