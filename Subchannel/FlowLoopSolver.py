"""
Solvers for flow loops.
"""

# NUMPY
import numpy as np

# PANDAS
import csv

# MOTHMAN
from .FluidRelation import FluidRelation
from Meshing.Meshing import *
from Fields.Fields import *
from Kernels.Kernels import *
from Solvers.Solvers import *
from Subchannel.Channel import *

# MPL
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.patches import Rectangle

# PICKLE
import pickle as pkl

class FlowLoopTransientSolver:
  def __init__(self, start_time: float, end_time: float, Nsteps: int,
               channels: list, pre_solve_methods: list):

    self.start_time = start_time
    self.end_time = end_time
    self.Nsteps = Nsteps
    self.timesteps = np.linsapce(self.start_time, self.end_time, self.Nsteps)

    for chan in channels:
      if not isinstance(chan, Channel):
        raise TypeError("channels must be a list of Channel() types.")
    self.channels = channels

  def solve_transient(self):

    # UPDATE OLD SOLUTIONS TO MOST RECENT SOLUTIONS
    for chan in self.channels:
      chan.update_old_to_most_recent()


    # SET T PREV
    t_prev = -1e321
    for t_new in self.timesteps:
      # DT
      this_dt = t_new - t_prev

      # PRINT TIME INFORMATION
      print("NOW SOLVING AT TIME =", t_new, "| this_dt =", this_dt)

      # PRESOLVE METHODS
      for method in self.pre_solve_methods:


