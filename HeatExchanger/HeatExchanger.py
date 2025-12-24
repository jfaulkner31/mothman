"""
Imports
"""
import numpy as np
from Kernels.Helpers import *

"""
HX object
"""

class HeatExchangerObject:
  def __init__(self, area: float, htc: float,
               hot_dict: dict, cold_dict: dict,
               temp_rel_tol: float):
    """
    Labelling
    ==========
    cd : hot side
    cl : cold side
    """
    self.Tcd_in
    self.Tcd_out
    self.Tcl_in
    self.Tcl_out

    self.T_cd_in_old
    self.T_cl_in_old
    self.T_cd_out_old
    self.T_cl_out_old

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

  def set_hot_bc(self, T_in: float, mdot: float, pressure: float):
    self.T_cd_in = T_in
    self.hot['mdot'] = mdot
    self.hot['pressure'] = pressure

  def set_cold_bc(self, T_in: float, mdot: float, pressure: float):
    self.T_cl_in = T_in
    self.cold['mdot'] = mdot
    self.cold['pressure'] = pressure

  def update_old_to_most_recent(self):
    """
    Called externally, updates old to most recent.
    Usually called when soln. is done.
    """
    self.T_cd_in_old = self.T_cd_in
    self.T_cl_in_old = self.T_cl_in
    self.T_cd_out_old = self.T_cd_out
    self.T_cl_out_old = self.T_cl_out

  def solve(self, _dt: float):
    """
    Solves heat exchanger

    Params
    ======
    _dt : float
      time step in seconds
    """
    lmtd_old = self._logmean()
    diff = 1e6
    while diff > self.rel_tol:
      T_cd_out = self._hot_diffeq(_dt=self._dt, info=self.hot)
      T_cl_out = self._cold_diffeq(_dt=self._dt, info=self.cold)
      self.T_cd_out = T_cd_out
      self.T_cl_out = T_cl_out
      lmtd_new = self._logmean()
      diff = np.abs(1.0 - lmtd_old/lmtd_new)
      lmtd_old = lmtd_new

  def _logmean(self):
    """
    Returns logmean temperature difference based on
    current temperature values.

    Params
    ======
    None
    """
    top = self.T_cd_in - self.T_cl_in - (self.T_cd_out - self.T_cl_out)
    bottom = np.log(self.T_cd_in - self.T_cl_in) - np.log(self.T_cd_out - self.T_cl_out)
    return top / bottom

  def _heat_source_term(self, info: dict):
    M =  info['mass']
    cp = info['cp']
    return 2.0 / M / cp * self.get_H() * self.A * self.htc * self._logmean()

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
    C1 = 2.0 * mdot * M
    C2 = self._heat_source_term(info=info)

    Tdiff = (self.T_cd_in - self.T_cd_in_old) * F

    T_cd_out = (-C1 * self.T_cd_in + C2 - Tdiff) / (F - C1)
    return T_cd_out

  def _cold_diffeq(self, _dt: float, info: dict) -> float:
    """
    Differential equation for the cold side
    """
    M =  info['mass']
    mdot = info['mdot']
    cp = info['cp']

    F = 1.0/_dt
    C1 = 2.0 * mdot * M
    C2 = self._heat_source_term(info=info)

    Tdiff = (self.T_cd_in - self.T_cd_in_old) * F

    T_cl_out = (-C1 * self.T_cd_in - C2 - Tdiff) / (F - C1)
    return T_cl_out



