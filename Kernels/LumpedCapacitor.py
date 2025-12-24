from Meshing.Meshing import *
from Fields.Fields import *
from Kernels.Kernels import *
import matplotlib.pyplot as plt
from copy import deepcopy
from Kernels.Linkers import FloatLinker, Linker


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



class Conductor:
  """
  Conductor base class for all conductors.

  Params
  ======
  A : float
    wall area for heat transfer
  wall_temp : float
    temperature of the wall for conjugate heat transfer

  Methods
  =======
  Methods for the default class should include any methods
  that must be called by the subchannel application -> e.g.
  all methods much have a get_wall_area() and get_wall_temp() method
  for retrieving the wall area and wall temperature respectively.

  solve() : method that performs the solution of the conductor
  get_wall_area() : returns wall area for heat transfer
  get_wall_temp() : returns wall temperature for heat transfer
  """

  def __init__(self):
    self.A = None
    self.wall_temp = None
    self.power = None
  def solve(self):
    pass
  def get_wall_area(self) -> float:
    return self.A
  def get_wall_temp(self) -> float:
    return self.wall_temp
  def get_max_temp(self) -> float:
    raise Exception("No default get_max_temp implemented!")
  def get_min_temp(self) -> float:
    raise Exception("No default get_min_temp implemented!")
  def get_mean_temp(self) -> float:
    raise Exception("No default get_min_temp implemented!")
  def get_integrated_power(self) -> float:
    """method for getting integrated power int(P''' dV) -> float"""
    raise Exception("No default method for getting integrated power!")
  def get_coupling_coeffs() -> tuple[np.ndarray, np.ndarray, float]:
    raise Exception("No default method for getting coupling coeffs!")
  def set_T(self, T: float | np.ndarray) -> None:
    raise Exception("No default method for setting T!")


class LumpedCapacitor(Conductor):
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
  def __init__(self, mass: float, power: float | Linker, h: float | Linker, A: float,
               C: float | LumpedMaterialProperty, thermal_cond: float | LumpedMaterialProperty,
               initial_T: float,
               T_bulk: float | Linker,
               epsilon: float):
    self.mass = mass # mass in kg
    self.power = power # Power in Watts
    self.h = h # htc in W/m2-K
    self.A = A # area in m2
    self.C = C # specific heat
    self.thermal_cond = thermal_cond # thermal cond
    self.T_bulk = T_bulk # bulk temperature - use a Linker or a float

    self.T = initial_T
    self.T_old = self.T

    self.epsilon = epsilon

  def set_T(self, T: float):
    """
    Sets the value of Temperature (K)
    """
    self.T = T

  def solve(self, _dt: float):
    """
    Solves

    T_next = k+1
    self.T = k
    T_old  = k-1
    """

    T_bulk = self.T_bulk

    # If it is a lumped material property.
    if isinstance(self.C, LumpedMaterialProperty):
      C = self.C.explicit_eval(self.T)
      iterate = not self.C.is_constant()
    else:
      C = self.C
      iterate = False

    if iterate:
      diff = 1e321
      T_next_prev_guess = 1e321
      while (diff > self.epsilon):
        # Double check tihs - specifically if self.T in the htc term should actually be
        # part of T_next -> use the forward scheme or do we use the old value of T?
        # or do we use T_old?
        # Getting my T's mixed up a lot here.
        T_next = self.T_old + _dt * (
          self.power / self.mass / C + self.h*self.A / self.mass / C * (T_bulk - self.T)
        )
        diff = np.abs(T_next - T_next_prev_guess)
        T_next_prev_guess = T_next
    else:
      # T_next = self.T_old + _dt * (
      #   self.power / self.mass / C + self.h*self.A / self.mass / C * (T_bulk - self.T)
      # )
      bottom_term = 1/_dt + self.h*self.A / self.mass / C
      top = self.power / self.mass / C + self.h*self.A / self.mass / C * T_bulk + self.T_old / _dt
      T_next = top / bottom_term

    # Save old and new T
    #self.T_old = copy.deepcopy(self.T)
    self.T = T_next


  def get_coupling_coeffs(self, _dt: float)->tuple[np.ndarray, np.ndarray, float]:
    """
    Basic function for doing a coupled
    solve with a channel object.

    Returns matrix H and vector B as objects

    Returns coupling coefficient as float
    where coupling coefficient is to be multiplied by T_bulk

    Does not consider the heat flux term: HAT_b / MCp

    H = 1/_deltaT + HA/MC -> matrix
    B = P/MC + T_old/_deltaT -> vector
    coeff = HA/MC -> float

    """

    # If it is a lumped material property.
    if isinstance(self.C, LumpedMaterialProperty):
      C = self.C.explicit_eval(self.T)
    else:
      C = self.C

    H = np.array([1.0/_dt + self.h * self.A / self.mass / C])
    B = np.array([self.power / self.mass / C + self.T_old / _dt])
    coeff = self.h * self.A / self.mass / C * -1.0 # to be multiplied by T_bulk
    return H, B, coeff



  def update_old_to_most_recent(self):
    self.T_old = copy.deepcopy(self.T)

  ### GETTERS
  def get_heat_flux(self):
    """
    Returns heat flux:

      q'' = h*A*(T_f - T_solid)
    """
    return self.A * self.h * (self.T_bulk - self.T)

  def get_wall_temp(self):
    return self.T
  def get_wall_area(self):
    return self.A
  def get_max_temp(self):
    return self.T
  def get_min_temp(self):
    return self.T
  def get_mean_temp(self):
    return self.T
  def get_integrated_power(self) -> float:
    return self.power

  ### SETTERS
  def set_power(self, power: float | Linker):
    self.power = power

  def set_htc(self, h: float | Linker):
    self.h = h

  def set_T_bulk(self, T_bulk: float | Linker):
    self.T_bulk = T_bulk

  def set_fluid_vals(self, power: float, htc: float, T_bulk: float):
    self.set_power(power=power)
    self.set_htc(htc=htc)
    self.set_T_bulk(T_bulk=T_bulk)

  ### OTHER
  def update_thermal_props(self):
    """
    Updates/sets thermal_cond or C values if they are LumpedMaterialProperty
    """
    pass

