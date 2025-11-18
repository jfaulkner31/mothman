"""
This file contains helper functions when two things cannot be connected very nicely.
"""

"""
Imports
"""
from Subchannel.Channel import *
from Kernels.LumpedCapacitor import *

"""
Helper functions for connecting conductors to
channel.conductor objects.

Note that channel objects hold conductors -
ChannelArray objects do not.
"""

def ChannelArrayConductionBuilder(channelarray: ChannelArray,
               conductor: Conductor) -> None:
  """
  A class for helping setup a repetitive
  conductor for an array of channels

  Params
  ======
  channelarray : ChannelArray
    The channel array we are making conductors for
  conductor : Conductor
    The reference conductor object we are using.

  """
  for ch in channelarray.channels:
    # Make a ch.conductors list for this channel.
    the_arr: list[Conductor] = np.ndarray([])
    for idx in ch.temp.T:
      # Make a new conductor for this node
      the_arr = np.append(the_arr, copy.deepcopy(conductor))
    # Append the
    ch.conductors = the_arr[1:]

def ChannelConductionBuilder(channel: Channel,
                             conductor: Conductor) -> None:
  """
  A class for helping setup a repetitive
  conductor for a channel.

  Params
  ======
  channel : Channel
    The channel we are making conductors for
  conductor : Conductor
    The reference conductor object we are using.
  """
  the_arr: list[Conductor] = np.ndarray([])
  for idx in channel.temp.T:
    # Make a new conductor for this node
    the_arr = np.append(the_arr, copy.deepcopy(conductor))
  # Append the
  channel.conductors = the_arr[1:]
