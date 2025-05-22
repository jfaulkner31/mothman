from Meshing.Meshing import *
from Fields.Fields import *
from Kernels.Kernels import *
import matplotlib.pyplot as plt
from copy import deepcopy

class LumpedMaterialProperty:
  """
  Class that serves as a way to evaluate thermal properties.
  """
  def __init__(self):
    pass

class LumpedCapacitor:
  """
  Serves as a lumped capacitor for heat transfer:
  dT/dt = P/m/C + h*A/m/C * (T_bulk - T)
  T = solid temperature (solving for)
  P = power (W)
  m = mass (kg)
  C = specific heat (J/kg/K)
  A = heat transfer area to fluid
  """
  def __init__(self, mass: float, power: float, h: float, A: float, L: float,
               C: float | LumpedMaterialProperty, thermal_cond: float | LumpedMaterialProperty,
               initial_T: float,
               T_bulk: float):
    self.mass = mass
    self.power = power
    self.h = h
    self.A = A
    self.C = C
    self.thermal_cond = thermal_cond
    self.L = L
    self.T_bulk = T_bulk

    self.T = initial_T
    self.T_old = self.T

  def solve(self, _dt: float):
    """
    Solves.
    """
    T_next = self.T + _dt * (
      self.P / self.mass / self.C + self.h*self.A / self.mass / self.C * (self.T_bulk - self.T)
    )
    self.T_old = copy.deepcopy(self.T)
    self.T = T_next
  def get_heat_flux(self):
    """
    Returns heat flux.
    """
    return self.A * self.h * (self.T_bulk - self.T)
  def set_power(self, power: float):
    self.power = power
  def set_htc(self, h: float):
    self.h = h
  def set_T_bulk(self, T_bulk: float):
    self.T_bulk = T_bulk
  def update_thermal_props(self):
    """
    Updates/sets thermal_cond or C values if they are LumpedMaterialProperty
    """
    pass



