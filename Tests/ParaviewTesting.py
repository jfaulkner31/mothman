from Paraview.Mesher import ParaviewMeshFromChannelArray as PMesh
from Subchannel import *
import pickle as pkl

with open("../Data/1152_channels.pkl", "rb") as f:
  obj = pkl.load(f)

ch = obj['chArr']

dx, dy = 0.0254 * 2.0/2**0.5, 0.0254 * 2.0/2**0.5
pmesh = PMesh(arr=ch, filename='output.vtu', dx=dx, dy=dy)
pmesh.write()
