import numpy as np
from Meshing.Meshing import Mesh_1D
import copy

# General field class
class Field:
  def __init__(self, name: str, initial_value: (np.ndarray, float), mesh: Mesh_1D):
    self.name = name
    self.mesh = mesh
    self.centroids = self.mesh.centroids

# General scalar field class (a field with a single scalar variable such as Temperature)
class ScalarField(Field):
  """
  Scalar variable for transport model - represented by a T
  """
  def __init__(self, name: str, initial_value: (np.ndarray, float), mesh: Mesh_1D):
    # super from Field class
    super().__init__(initial_value=initial_value, name=name, mesh=mesh,)

    # Now to make initial values for this field based on if array or float
    if isinstance(initial_value, np.ndarray):
      if len(initial_value) == self.mesh.nz:
        self.initial_value = initial_value
      else:
        raise Exception("nz on variable input mesh must be the same as the length of the initial condition vector!")

    elif isinstance(initial_value, float):
      self.initial_value = np.ones(self.mesh.nz)*initial_value

    else:
      raise Exception("Initial value must be of np.ndarray or float.")

    # Assigns T as initial values based on deepcopy
    self.T = copy.deepcopy(self.initial_value)

    # More stuff..

class FaceField:
  """
    A field where values are stored on the faces.
  """
  def __init__(self, name: str, initial_value: (np.ndarray, float), mesh: Mesh_1D):
    self.name = name
    self.mesh = mesh
    self.centroids = self.mesh.centroids
    self.faceLocs = self.mesh.nodeCoords

    # Now to make initial values for this field based on if array or float
    if isinstance(initial_value, np.ndarray):
      if len(initial_value) == self.mesh.nz+1: # +1 because this is the face values.
        self.initial_value = initial_value
      else:
        raise Exception("nz on variable input mesh must be the same as the length of the initial condition vector!")

    elif isinstance(initial_value, float):
      self.initial_value = np.ones(self.mesh.nz+1)*initial_value

    else:
      raise Exception("Initial value must be of np.ndarray or float.")

    self.T = copy.deepcopy(self.initial_value)

  def get_upper(self, cid: int):
    """
      Takes in a cid for a cell and returns the upper face value.
    """
    return self.T[cid+1]
  def get_lower(self, cid: int):
    """
    Takes in a cid for a cell and returns the lower face value
    """
    return self.T[cid]
  def set_lower(self, cid: int, value: float):
    self.T[cid] = value
  def set_upper(self, cid: int, value: float):
    self.T[cid+1] = value
