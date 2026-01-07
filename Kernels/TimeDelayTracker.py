import numpy as np

"""
A basic class for handling time delays.
"""

class TimeDelayTracker:
  """
  Basic time delay tracker. Does no physics.
  Models V_out(t) = V(t - delay)
  where V is the tracked value.
  Linear interpolation is used for values not directly tracked.
  """
  def __init__(self, delay: float):
    self.times = np.array([])
    self.values = np.array([])
    self.delay = delay

  def add_pair(self, time: float, value: float):
    # First find if time is in times -> if it is, overwrite based on pos
    if time in self.times:
      pos = np.where(self.times == time)[0][0]
      self.values[pos] = value
    else:
      self.times = np.append(self.times, time)
      self.values =  np.append(self.values, value)

  def get_value(self, time: float) -> float:
    """
    Gets value at specified time. Linear interpolation is used if time is not
    directly tracked.
    """

    # If time is in times, return corresponding value
    if time in self.times:
      pos = np.where(self.times == time)[0][0]
      return self.values[pos]
    else:
      # Otherwise do linear interpolation
      if time < np.min(self.times) or time > np.max(self.times):
        raise Exception("Requested time is outside of tracked range")
      else:
        pos_right = np.where(self.times > time)[0][0]
        pos_left = pos_right - 1
        t_left = self.times[pos_left]
        t_right = self.times[pos_right]
        v_left = self.values[pos_left]
        v_right = self.values[pos_right]
        # Linear interpolation
        value = v_left + (v_right - v_left) * (time - t_left) / (t_right - t_left)
        return value

  def get_delayed_value(self, time: float) -> float:
    """
    Gets value at specified time minus delay

    Params
    ======
    time : float
      current time in seconds
    Returns
    =======
    value : float
      V_out(t) = V(t - delay)
    """
    delayed_time = time - self.delay
    return self.get_value(delayed_time)


  def clear_before(self, time: float):
    """
    Clears all tracked values before specified time
    """
    mask = self.times >= time
    self.times = self.times[mask]
    self.values = self.values[mask]

  def clear_before_minus_delay(self, time: float):
    """
    Clears all tracked values before (specified time - delay)
    """
    delayed_time = time - self.delay
    self.clear_before(delayed_time)

class ExponentialDelayTracker(TimeDelayTracker):
  """
  Exponential decay time delay tracker

  Params
  ======
  delay : float
    time delay amount in seconds
  decay : float
    exponential decay constant (1/seconds)

  Methods
  =======
  add_pair(time: float, value: float)
    adds a time, value pair
  get_value(time: float) -> float
    gets value at specified time that is stored in the dict
  get_delayed_value(time: float) -> float
    gets value at specified time minus delay, with exponential decay applied
  clear_before(time: float)
    clears all tracked values before specified time

  """
  def __init__(self, delay: float, decay: float):
    super().__init__(delay=delay)
    self.decay = decay

  def get_value(self, time: float) -> float:
    return super().get_value(time)

  def add_pair(self, time: float, value: float):
    return super().add_pair(time, value)

  def clear_before(self, time):
    return super().clear_before(time)

  def clear_before_minus_delay(self, time):
    return super().clear_before_minus_delay(time)

  def get_delayed_value(self, time):
    """
    Input the current time, returns the value at (time - delay)
    with exponential decay applied.

    Params
    ======
    time : float
      current time in seconds
    Returns
    =======
    value : float
      value at (time - delay) with exponential decay applied
      V_out(t) = V(t - delay) * exp(-decay * tau)
    """
    return self.get_value(time) * np.exp(-self.decay * self.delay)
