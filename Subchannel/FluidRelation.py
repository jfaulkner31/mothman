

class FluidRelation:
  def __init__(self, cp: float, mu: float, k: float, rho_0: float, drho_dT: float):
    self.cp = cp
    self.mu = mu
    self.k = k
    self.rho_0 = rho_0 # example 2715.13 - 0.513 * T[K] ---> rho_eos = [2715.13, -0.513] -> rho_0 = 2715.13 ; drho_dT = -0.513
    self.drho_dT = drho_dT
  def props_from_P_H(self, P: float, enthalpy: float, prop: str):
    # First get temperature from pressure and enthalpy
    _T = enthalpy / self.cp # ok
    if prop == 'T':
      return _T
    elif prop == 'rho':
      return self.rho_0 + self.drho_dT * _T
    else:
      raise Exception("Requested property not available")

  def props_from_P_T(self, P:float, T:float, prop: str):
    if prop == 'h':
      return self.cp*T
    elif prop == 'rho':
      return self.rho_0 + self.drho_dT * T
    else:
      raise Exception("Desired prop not found!")
  def get_mu(self):
    return self.mu
  def get_k(self):
    return self.k
