from Meshing.Meshing import *
from Fields.Fields import *
from Kernels.Kernels import *
import matplotlib.pyplot as plt
from copy import deepcopy

class LumpedMaterialProperty:
  """
  Class that serves as a way to evaluate thermal properties.

  name: material property name

  function_coeffs: assumes a polynomial of the form:
  fc[0] + fc[1]*x + fc[2]*x^2 + fc[3]*x^3 ...
  """
  def __init__(self, name: str, function_coeffs: list | np.ndarray):
    self.fc = function_coeffs
    self.name = name

  def explicit_eval(self, x: float | np.ndarray):
    """
    Input x and returns material property value based on function coeffs.
    """
    if isinstance(x, float):
      power = 0.0
      result = 0.0
      for this in self.fc:
        result += this * x**power
        power += 1.0
      return result

    elif isinstance(x, np.ndarray):
      result = np.zeros(len(x))
      power = 0.0
      for this in self.fc:
        result += this * x**power
        power += 1.0
      return result

    else:
      raise TypeError("x input is wrong type!")
  def plot_prop(self, xMin: float, xMax: float):
    xRange = np.linspace(xMin,xMax, 1000)
    y = self.explicit_eval(x=xRange)
    ax = plt.figure()
    plt.plot(xRange, y, 'k-', linewidth=2)
    plt.plot(xRange[0::99], y[0::99], 'rs', markerfacecolor='w', markersize=6)
    plt.xlabel('x', fontsize=15)
    plt.ylabel(self.name, fontsize=15)
    plt.grid()
    return ax


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
               T_bulk: float,
               epsilon: float):
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

    self.epsilon = epsilon

  def solve(self, _dt: float):
    """
    Solves.
    """
    if isinstance(self.C, LumpedMaterialProperty):
      C = self.C.explicit_eval(self.T)
      iterate = True
    else:
      C = self.C
      iterate = False

    if iterate:
      diff = 1e321
      T_next_prev_guess = 1e321
      while (diff > self.epsilon):
        T_next = self.T + _dt * (
          self.P / self.mass / C + self.h*self.A / self.mass / C * (self.T_bulk - self.T)
        )
        diff = np.abs(T_next - T_next_prev_guess)
        T_next_prev_guess = T_next
    else:
      T_next = self.T + _dt * (
        self.P / self.mass / C + self.h*self.A / self.mass / C * (self.T_bulk - self.T)
      )

    # Save old and new T
    self.T_old = copy.deepcopy(self.T)
    self.T = T_next

  def get_heat_flux(self):
    """
    Returns heat flux:

      q'' = h*A*(T_f - T_solid)
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



