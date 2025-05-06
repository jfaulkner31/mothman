import numpy as np
from Meshing.Meshing import *
from Fields.Fields import *
from Kernels.Kernels import *
from Solvers.Solvers import *
from Subchannel.FluidRelation import FluidRelation
from Subchannel.Channel import Channel
from Subchannel.Channel import ChannelInterface
from Subchannel.Channel import ChannelArray
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
  integral = sum(C * lam * channel.vol_vec * weight)

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

    # Get physical delayed neutron fraction
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

def compute_beff_multichannel(channels: list, names: list, weights: list):
  beta_num_list = []
  num_int = 0.0
  beta_sum = 0.0

  # Integrate Numerator
  for idx, this_ch in enumerate(channels):
    for name_idx, this_name in enumerate(names):
      # CHANNELS
      this_int = compute_beff_numerator(channel=this_ch, name=this_name, weight=weights[idx])
      num_int += this_int
      # STORE RESULTS
      try:
        beta_num_list[name_idx] += this_int
      except:
        beta_num_list.append(this_int)

  # Get physical DNP's for each channel
  beta_physical_by_channel = []
  for this_ch in channels:
    beta_phys_sum = 0.0
    for name in names:
      beta_phys_sum += this_ch.tracer_kernels[name][2].beta
    beta_physical_by_channel.append(beta_phys_sum)

  # Integrate denominator
  prompt_int = 0.0
  for idx, this_ch in enumerate(channels):
    beta_sum_this_ch = beta_physical_by_channel[idx]

    # Get weight as a field
    nz = len(this_ch.mesh.cidList)
    if isinstance(weights[idx],float):
      this_wt = np.ones(nz)*weights[idx]
    elif isinstance(weights[idx], ScalarField):
      this_wt = copy.deepcopy(weights[idx].T)
    else:
      this_wt = copy.deepcopy(weights[idx])

    # denominator math now
    prompt_int += sum(this_wt * this_wt * this_ch.vol_vec) # * (1-beta_sum_this_ch)

  # Get real beff
  beff = num_int / (num_int + prompt_int)

  # Beff by dnp:
  beta_eff_by_dnp = [this / (num_int + prompt_int) for this in beta_num_list]

  return beff, beta_eff_by_dnp
