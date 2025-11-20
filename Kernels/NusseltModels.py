"""
Nusselt number models for
getting a heat transfer coefficient.
"""

class NusseltModel:
  def __init__(self):
    pass
  def get_Nu(self) -> float:
    raise Exception("No default nusselt model implemented!")

class BasicNusseltModel(NusseltModel):
  def __init__(self):
    super().__init__()
  def get_Nu(self) -> float:
    return 48.0 / 11.0
