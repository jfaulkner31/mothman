import numpy as np
from Meshing.Meshing import *
from Fields.Fields import *
from Kernels.Kernels import *
from Solvers.Solvers import *
from Subchannel.FluidRelation import FluidRelation
from Subchannel.Channel import Channel
from Subchannel.Channel import ChannelInterface
import copy

def compute_beff_numerator(channel: Channel, name: str, weight: np.ndarray | float | ScalarField):
  # First get some values:
  C = channel.tracers[name].T
  lam = channel.tracer_kernels[name][1].lam
  beta = channel.tracer_kernels[name][2].beta

  # Now convert weight function to something usable:
  nz = len(channel.mesh.cidList)
  if isinstance(weight,float):
    weight = np.ones(nz)*weight
  elif isinstance(weight, ScalarField):
    weight = copy.deepcopy(weight.T)
  else:
    weight = weight

  # Now integrate:
  integral = 0.0
  for cid in channel.mesh.cidList:
    vol = channel.mesh.cells[cid].vol
    integral += C[cid] * lam * vol * weight[cid]

  return integral

def compute_beff(channel, names: list, weight: np.ndarray | float | ScalarField):
  beff_num_list = []
  num_int = 0.0
  beta_sum = 0.0

  # Integrate numerator
  for name in names:
    this_int = compute_beff_numerator(channel=channel, name=name, weight=weight)
    beff_num_list.append(this_int)
    num_int += this_int

    # Compute physicl delayed neutron fraction
    beta_sum += channel.tracer_kernels[name][2].beta

  # Now convert weight function to something usable:
  nz = len(channel.mesh.cidList)
  if isinstance(weight,float):
    weight = np.ones(nz)*weight
  elif isinstance(weight, ScalarField):
    weight = copy.deepcopy(weight.T)
  else:
    weight = weight

  # Integrate denominator Part 2 (Prompt term)
  prompt_int = 0.0
  for cid in channel.mesh.cidList:
    vol = channel.mesh.cells[cid].vol
    prompt_int += vol * weight[cid] * weight[cid] * (1-beta_sum)

  # Now add:
  denom = prompt_int + num_int
  delayed_int = sum(beff_num_list)
  beff = delayed_int / denom
  beta_i = np.array(beff_num_list) / denom

  return beff, beta_i



