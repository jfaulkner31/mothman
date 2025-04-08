import numpy as np
from Meshing.Meshing import Mesh_1D
import copy

class ScalarField:
  """
  Scalar variable for transport model - represented by a T
  """
  def __init__(self, name: str, initial_value: (np.ndarray, float), mesh: Mesh_1D):
    self.name = name
    self.mesh = mesh
    if isinstance(initial_value, np.array):
      if len(initial_value) == self.mesh.nz:
        self.initial_value = initial_value
      else:
        raise Exception("nz on variable input mesh must be the same as the length of the initial condition vector!")
    elif isinstance(initial_value, float):
      self.initial_value = np.ones(self.mesh.nz)*initial_value
    else:
      raise Exception("Initial value must be of np.ndarray or float.")
    self.T = copy.deepcopy(self.initial_value)
