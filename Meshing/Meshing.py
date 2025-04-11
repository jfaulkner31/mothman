import numpy as np

class Cell:
  def __init__(self, dz: float, upperZ: float, lowerZ: float, upperType: str, lowerType: str, upperArea: float, lowerArea: float, cid: int, upperFid: int, lowerFid: int):
    self.cid = cid # cell id
    self.upperZ = upperZ # upper face position
    self.lowerZ = lowerZ # lower face position
    self.dz = self.upperZ - self.lowerZ # size of cell aka the cell delta value
    self.centroid = (self.upperZ + self.lowerZ) / 2.0
    self.upperType = upperType # type of boundary - f for internal face and b for boundary
    self.lowerType = lowerType
    self.upperArea = upperArea # upper face area
    self.lowerArea = lowerArea # lower face area
    self.dFaceUpper = self.upperZ - self.centroid # absolute distance to upper FACE
    self.dFaceLower = self.centroid - self.lowerZ # absolute distance to lower FACE

    self.dCellUpper = None # absolute distance to lower or upper face
    self.dCellLower = None

    self.lowerNeighborCid = None # stays at None if no lower neighbor - otherwise is lower neighbor cell id
    self.upperNeighborCid = None

    self.lowerNeighborDist = None # stays at None if no dist - otherwise is lower neighbor distance to centroid
    self.upperNeightborDist = None

    self.upperNeighborGdiff = None
    self.lowerNeighborGdiff = None

    self.upperFid=upperFid
    self.lowerFid=lowerFid

class Mesh_1D:
  """
  Creates a one dimensional mesh object.
  """
  def __init__(self, nodeCoords: np.ndarray, faceAreas: np.ndarray):

    # proofreading / exceptions
    self.nodeCoords = np.unique(np.sort(nodeCoords)) # nodeCoords = locations of zFaces starting at boundary
    if not np.array_equal(self.nodeCoords, nodeCoords):
      raise Exception("Node coords are not sorted or has repeating values in it!")
    if len(faceAreas) != len(self.nodeCoords):
      raise Exception("Face areas are not the same length as node coords input!")
    self.faceAreas = faceAreas
    self.cells = {} # dict of cells  - key is cell id and value is the cell object
    self.nz = len(self.nodeCoords)-1 # number
    self.cidList = []
    self.centroids = []

    # now assign cells and stuff
    cid = int(0)
    fid = int(0)
    for idx, node in enumerate(self.nodeCoords):
      if idx == 0:
        continue
      else:
        upperType = 'f'
        lowerType = 'f'
        if idx == len(self.nodeCoords)-1: # last so upper face is a boundary
          upperType = 'b'
        if idx == 1: # lower face is a boundary
          lowerType = 'b'

        # assign cell and its data
        self.cells[cid] = Cell(cid=cid,
                               dz=self.nodeCoords[idx] - self.nodeCoords[idx-1],
                               upperZ=self.nodeCoords[idx],
                               lowerZ=self.nodeCoords[idx-1],
                               upperType=upperType,
                               lowerType=lowerType,
                               upperArea=self.faceAreas[idx],
                               lowerArea=self.faceAreas[idx-1],
                               upperFid=fid,
                               lowerFid=fid-1)
        self.cidList.append(cid)
        self.fidList.append(fid)


        fid += 1 # progress fid
        cid += 1 # progress cid

    # now we iterate across all cells and define distance to neighbors and neighbor ids
    # ids and distances will be None if the neighbor does not exist
    for cid in self.cells.keys():
      this_cell = self.cells[cid]
      # lower cell try
      try:
        this_cell.lowerNeighborCid = self.cells[this_cell.cid - 1].cid
        this_cell.lowerNeighborDist = abs(this_cell.centroid - self.cells[this_cell.cid - 1].centroid)
        this_cell.lowerNeighborGdiff = this_cell.lowerArea / this_cell.lowerNeighborDist
      except:
        pass

      # upper cell try
      try:
        this_cell.upperNeighborCid = self.cells[this_cell.cid + 1].cid
        this_cell.upperNeighborDist = abs(this_cell.centroid - self.cells[this_cell.cid + 1].centroid)
        this_cell.upperNeighborGdiff = this_cell.upperArea/this_cell.upperNeighborDist
      except:
        pass

    # finally make a z vector
    for cid in self.cells.keys():
      self.centroids.append(self.cells[cid].centroid)




