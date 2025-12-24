"""
This file contains helper functions and classes.

Examples:
  Connecting a conductor to all channels in a ChannelArray.
  Connect a conductor to a Channel.

"""

"""
Imports
"""
from Subchannel.Channel import *
from Kernels.Linkers import Linker, FloatLinker

"""
Helper function for setting integrated power of each.
"""
def ChannelSetTotalPower(ch: Channel, cond_power: float, channel_power: float):
  """
  Scales total channel power and the conductor powers so that:

  int(P''' dV)_conductor = cond_power
  int(P''' dV)_channel   = channel_power

  Params
  ======
  channel : Channel
    The channel
  cond_power : float
    conductor integrated power in Watts
  channel_power : float
    channel power we want to scale to.
  """

  # Get original powers
  orig_ch_power =   ch.get_integrated_power()
  orig_cond_power = ch.get_integrated_conductor_power()

  cond_ratio = cond_power / orig_cond_power
  fluid_ratio = channel_power / orig_ch_power

  # First scale the conductor powers
  for cond in ch.conductors:
    if isinstance(cond.power, FloatLinker):
      # Is a linker so we must keep it as a linker.
      # Multiply linker multiplier by ratio of OLD to NEW
      if cond.power.multiplier == 0:
        raise Exception("Multiplier is == 0 -> exception!")
      cond.power.multiplier *= cond_ratio / fluid_ratio

    elif isinstance(cond.power, float):
      # Is a float
      cond.power *= cond_ratio / fluid_ratio

    elif isinstance(cond.power, list):
      raise Exception("List option not yet implemented for power scaling!")

  # Next scale channel power
  ch.heat_source.T *= fluid_ratio

  # Get stuff new
  new_cond_power = ch.get_integrated_conductor_power()
  new_ch_power   = ch.get_integrated_power()
  total_power    = new_ch_power + new_cond_power
  # Printing
  print(f'Old conductor power = {orig_cond_power}')
  print(f'New conductor power = {new_cond_power}')
  print(f'Old channel power = {orig_ch_power}')
  print(f'New channel power = {new_ch_power}')
  print(f'Total new power = {total_power}')
  print(f'Conductor power fraction = {new_cond_power/total_power}')
  print(f'Channel power fraction = {new_ch_power/total_power}')

def ChannelArraySetTotalPower(arr: ChannelArray, cond_power: float, channel_power: float) -> None:
  """
  Scales total channel power and the conductor powers so that:

  int(P''' dV)_conductor = cond_power
  int(P''' dV)_channel   = channel_power

  Params
  ======
  channel : ChannelArray
    The channel array
  cond_power : float
    conductor integrated power in Watts
  channel_power : float
    channel power we want to scale to.
  """

  # Get original powers
  orig_ch_power =   arr.get_integrated_power()
  orig_cond_power = arr.get_integrateds_conductor_power()

  cond_ratio = cond_power / orig_cond_power
  fluid_ratio = channel_power / orig_ch_power

  # First scale the conductor powers
  for ch in arr.channels:
    for cond in ch.conductors:
      if isinstance(cond.power, FloatLinker):
        # Is a linker so we must keep it as a linker.
        # Multiply linker multiplier by ratio of OLD to NEW
        if cond.power.multiplier == 0:
          raise Exception("Multiplier is == 0 -> exception!")
        cond.power.multiplier *= cond_ratio / fluid_ratio

      elif isinstance(cond.power, float):
        # Is a float
        cond.power *= cond_ratio / fluid_ratio

      elif isinstance(cond.power, list):
        raise Exception("List option not yet implemented for power scaling!")

  # Next scale channel power
  for ch in arr.channels:
    ch.heat_source.T *= fluid_ratio

  # Get stuff new
  new_cond_power = arr.get_integrated_conductor_power()
  new_ch_power   = arr.get_integrated_power()
  total_power    = new_ch_power + new_cond_power
  # Printing
  print(f'Old conductor power = {orig_cond_power}')
  print(f'New conductor power = {new_cond_power}')
  print(f'Old channel power = {orig_ch_power}')
  print(f'New channel power = {new_ch_power}')
  print(f'Total new power = {total_power}')
  print(f'Conductor power fraction = {new_cond_power/total_power}')
  print(f'Channel power fraction = {new_ch_power/total_power}')




"""
Helper functions for connecting conductors to
channel.conductor objects.

Note that channel objects hold conductors -
ChannelArray objects do not.
"""

def ChannelArrayConductionBuilder(channelarray: ChannelArray,
               conductor: Conductor, ch_idx_DNI) -> None:
  """
  A class for helping setup a repetitive
  conductor for an array of channels.

  Copies the channels

  Params
  ======
  channelarray : ChannelArray
    The channel array we are making conductors for
  conductor : Conductor
    The reference conductor object we are using.
  ch_idx_DNI : list[int]
    The list of channel indices that we dont add the conductors to.

  """
  for idx, ch in enumerate(channelarray.channels):
    # Skip the channels we do not want to add conductors to.
    if idx in ch_idx_DNI:
      continue

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

"""
Helper functions for quickly linking
channel powers and their conductor powers.

Makes it so that when the channel powers are set/changed,
the powers in the graphite are also automatically set/changed.

AKA links Lumped capacitor power to channel heat source power.
"""

def LinkConductorPowersToChannel(channel: Channel, multiplier: float):
  """
  Links the power values of the conductors to the powers in the channels
  via a FloatLinker with a ratio chosen.
  """
  for idx, cond in enumerate(channel.conductors):
    cond.set_power(FloatLinker(idx=idx, multiplier=multiplier, obj=channel.heat_source.T))

def LinkConductorPowersToChannelArray(channelarray: ChannelArray, multiplier: float,
                                 ch_idx_DNI : list | np.ndarray) -> None:
  """
  Links the power values of the conductors to the powers in the channels
  via a FloatLinker with a ratio chosen.

  Skips channels in the ch_idx_DNI list.
  """
  for ch_idx, channel in enumerate(channelarray.channels):
    # Skip the channels we do not want to add conductors to.
    if ch_idx in ch_idx_DNI:
      continue
    # Now we call the above function to link the powers
    LinkConductorPowersToChannel(channel=channel, multiplier=multiplier)

"""
Linking functions for channel htc to conductor htc's
"""
def LinkConductorHTCToChannel(channel: Channel, multiplier: float):
  """
  Links the htc values of the conductors to the htc values in the channels
  via an optional multiplier.
  """
  for idx, cond in enumerate(channel.conductors):
    cond.set_htc(FloatLinker(idx=idx, multiplier=multiplier, obj=channel.htc.T))

def LinkConductorHTCToChannelArray(channelarray: ChannelArray, multiplier: float,
                                 ch_idx_DNI : list | np.ndarray) -> None:
  """
  Links the htc values of the conductors to the htc values in the channel array channels
  via an optional multiplier.

  Skips channels in the ch_idx_DNI list.
  """
  for ch_idx, channel in enumerate(channelarray.channels):
    # Skip the channels we do not want to add conductors to.
    if ch_idx in ch_idx_DNI:
      continue
    # Now we call the above function to link the powers
    LinkConductorHTCToChannel(channel=channel, multiplier=multiplier)

"""
Linking functions for channel T_bulk to conductor T_bulk's
"""
def LinkConductorTbulkToChannel(channel: Channel, multiplier: float):
  """
  Links the T_bulk values of the conductors to the T_bulk values in the channels
  via an optional multiplier.
  """
  for idx, cond in enumerate(channel.conductors):
    cond.set_T_bulk(FloatLinker(idx=idx, multiplier=multiplier, obj=channel.temp.T))

def LinkConductorTbulkToChannelArray(channelarray: ChannelArray, multiplier: float,
                                     ch_idx_DNI : list | np.ndarray) -> None:
  """
  Links T_bulk values of the conductors to the T_bulk values in the channel array channels
  """
  for ch_idx, channel in enumerate(channelarray.channels):
    # Skip the channels we do not want to add conductors to.
    if ch_idx in ch_idx_DNI:
      continue
    # Now we call the above function to link the powers
    LinkConductorTbulkToChannel(channel=channel, multiplier=multiplier)
