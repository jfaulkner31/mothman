from Meshing.Meshing import *
from Fields.Fields import *
from Kernels.Kernels import *
import matplotlib.pyplot as plt
from copy import deepcopy
from Subchannel.Channel import *

class LumpedMaterialProperty:
  """
  Class that serves as a way to evaluate thermal properties.

  This class does not ever hold field or solution data, it should
  in theory be usable across as many materials / objects as possible.

  name: material property name

  function_coeffs: assumes a polynomial of the form:
  fc[0] + fc[1]*x + fc[2]*x^2 + fc[3]*x^3 ...
  """
  def __init__(self, name: str, function_coeffs: list | np.ndarray):
    # Basic stuff
    self.fc = function_coeffs # coeffs for polynomial
    self.name = name # name of variable
    if (len(self.fc) == 1):
      self._is_constant = True
    else:
      self._is_constant = False

    # Setup derivatives
    self.dc = np.zeros(len(self.fc)) # derivative coeffs [0 + 1*fc[1] + 2*fc[2]*x + 3*fc[3]*x^2 + ...]
    self._get_derivatives()

  def _get_derivatives(self) -> None:
    # Get derivatives
    for idx, term in enumerate(self.fc):
      self.dc = idx * self.fc[0]


  def explicit_eval(self, x: float | np.ndarray):
    """
    Input x and returns material property value based on function coeffs.
    """
    if isinstance(x, float) | isinstance(x, int):
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

  def plot_prop(self, xMin: float, xMax: float) -> None:
    xRange = np.linspace(xMin,xMax, 1000)
    y = self.explicit_eval(x=xRange)
    ax = plt.figure()
    plt.plot(xRange, y, 'k-', linewidth=2)
    plt.plot(xRange[0::99], y[0::99], 'rs', markerfacecolor='w', markersize=6)
    plt.xlabel('x', fontsize=15)
    plt.ylabel(self.name, fontsize=15)
    plt.grid()
    # return ax

  def is_constant(self) -> bool:
    return self._is_constant

class BulkObjectHelper:
  """
  Object representing some bulk fluid with some temeperature value.
  Not meant to be used as a class - make subclasses of this guy instead
  where each has different implementations for grabbing the current
  temperature value.

  Really just a class added to LumpedCapacitory
  to conveniently pass information from channel objects.
  """
  def __init__(self):
    pass
  def get_T():
    pass
  def get_T_old():
    pass
  def get_htc(self) -> None:
    raise Exception("method not yet implemented!!")

class ChannelTempLinker(BulkObjectHelper):
  """
  Subchannel node bulk object that will retrieve
  temperature values from a Channel
  """
  def __init__(self, node_index: int, ch: Channel):
    self.ch = ch
    self.idx = node_index

  def get_T(self) -> float:
    return self.ch.temp.T
  def get_T_old(self) -> float:
    return self.ch.temp.T_old
  def get_htc(self) -> None:
    raise Exception("method not yet implemented!!")


class LumpedCapacitor:
  """
  Serves as a lumped capacitor for heat transfer:

  Params
  ======
  mass : float
    mass in kg
  power : float
    power in W
  A : float
    heat transfer area to fluid in m2
  L : length
    unit length in m
  C : float | LumpedMaterialProperty
    specific heat in J/kg/K
  thermal_cond : float | LumpedMaterialProperty
    thermal conductivity of solid
  initial_T : float
    initial condition in K
  T_bulk : float
    bulk temperature of fluid
  epsilon : float
    convergence tolerance in K

  Physics
  =======
  dT/dt = P/m/C + h*A/m/C * (T_bulk - T)
  """
  def __init__(self, mass: float, power: float, h: float, A: float, L: float,
               C: float | LumpedMaterialProperty, thermal_cond: float | LumpedMaterialProperty,
               initial_T: float,
               T_bulk: float | BulkObjectHelper,
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
    Solves

    T_next = k+1
    self.T = k
    T_old  = k-1
    """
    T_bulk = None
    # Determine how we are getting the bulk temperature
    if isinstance(self.T_bulk, float) | isinstance(self.T_bulk, int):
      T_bulk = float(self.T_bulk)
    elif isinstance(self.T_bulk, BulkObjectHelper):
      T_bulk = self.T_bulk.get_T()

    # If it is a lumped material property.
    if isinstance(self.C, LumpedMaterialProperty):
      C = self.C.explicit_eval(self.T)
      iterate = not self.C.is_constant()
    else:
      C = self.C
      iterate = False
    iterate = True

    if iterate:
      diff = 1e321
      T_next_prev_guess = 1e321
      while (diff > self.epsilon):
        T_next = self.T + _dt * (
          self.power / self.mass / C + self.h*self.A / self.mass / C * (T_bulk - self.T)
        )
        diff = np.abs(T_next - T_next_prev_guess)
        T_next_prev_guess = T_next
    else:
      T_next = self.T + _dt * (
        self.power / self.mass / C + self.h*self.A / self.mass / C * (T_bulk - self.T)
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

  def set_power(self, power: float | BulkObjectHelper):
    self.power = power

  def set_htc(self, h: float | BulkObjectHelper):
    self.h = h

  def set_T_bulk(self, T_bulk: float | BulkObjectHelper):
    self.T_bulk = T_bulk

  def update_thermal_props(self):
    """
    Updates/sets thermal_cond or C values if they are LumpedMaterialProperty
    """
    pass



