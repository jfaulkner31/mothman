"""
Imports
"""
import numpy as np


"""
Helper classes for connecting values in simulation immutable objects to other values.

e.g. power linker could link values in channel heat source to power in lumped capacitors.

"""

class Linker:
  """
  Base class for linking values.
  """
  def __init__(self):
    pass
  def value(self) -> float:
    raise Exception("No default linker implemented!")


class FloatLinker(Linker):
  """
  Links a float in one place (via indexing) in order to return a float
  value when doing simple match operations. Multiplier is used to scale the linked float.

  Used for simple retreival of float values
  """
  def __init__(self, idx: int, multiplier: float, obj: np.ndarray):
    """
    Params
    ======
    idx : int
      The index of the float we are linking to
    ratio : float
      The ratio we are multiplying the linked float by
    """
    self.idx = idx
    self.multiplier = multiplier
    self.obj = obj # object reference that holds the float value

  """External function to get the value"""
  def value(self) -> float:
    return self._get_value()

  """Internal function to get the value"""
  def _get_value(self) -> float:
    return self.obj[self.idx] * self.multiplier

  """Operator overloads for basic math operations"""
  def __add__(self, other: float) -> float:
    return self._get_value() + other
  def __radd__(self, other: float) -> float:
    return self._get_value() + other
  def __sub__(self, other: float) -> float:
    return self._get_value() - other
  def __rsub__(self, other: float) -> float:
    return other - self._get_value()
  def __mul__(self, other: float) -> float:
    return self._get_value() * other
  def __rmul__(self, other: float) -> float:
    return self._get_value() * other
  def __truediv__(self, other: float) -> float:
    return self._get_value() / other
  def __rtruediv__(self, other: float) -> float:
    return other / self._get_value()
