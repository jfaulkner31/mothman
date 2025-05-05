import pandas as pd
import numpy as np

def load_csv(filepath):
  """Loads the CSV file into a pandas DataFrame."""
  return pd.read_csv(filepath)

def csv_interpolator(csv_df, x_value: float, x_label: str, y_label: str):
  """
  Returns a y value based on x-y values from a csv file.
  NOTE: Assumes that the x-y pairs are presorted if going by time dependent function (e.g. x = time; y = f(t))
  """
  x = csv_df[x_label].values
  y = csv_df[y_label].values
  if x_value > x[-1]:
    return y[-1]
  elif x_value <= x[0]:
    return y[0]
  else:
    return np.interp(x_value, x, y)
