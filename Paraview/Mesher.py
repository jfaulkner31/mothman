from Subchannel.Channel import Channel, ChannelArray
import xml.etree.ElementTree as ET
import numpy as np
import math

"""
Generic default class for paraview
"""
class ParaviewMesh:
  def __init__(self, filename: str):
    self.filename = filename
  def write(self, time: float):
    pass


"""
Paraview mesh from channelarray
"""
class ParaviewMeshFromChannelArray(ParaviewMesh):
  def __init__(self, arr: ChannelArray, filename: str, dx: float, dy: float, writeConductors: bool = False):
    self.arr = arr
    self.name = filename
    self.dx = dx
    self.dy = dy
    self.writeConductors = writeConductors

  def write(self):
    """
    arr : ChannelArray
      channel array object we want to dump
    dx : float
      length of rect. prism each channel has (x)
    dy : float
      length of rect. prism each channel has (y)
    """
    arr = self.arr
    dx = self.dx
    dy = self.dy
    writeConductors = self.writeConductors

    points = []
    cells = []
    cell_types = []
    offsets = []
    point_index = 0
    field_array = []
    fields_dict = {}

    if writeConductors:
      fields_dict['T_c_mean'] = []

    # Iterate over channels
    for idx, ch in enumerate(arr.channels):
      x_offset = ch.xCoord
      y_offset = ch.yCoord
      if (x_offset is None) | (y_offset is None):
        print(f"Skipping channel number {idx} becuase x or y coord is None")
        continue
      x_offset /= 100.0
      y_offset /= 100.0

      # Setup the fields dict if ont set up
      the_fields = ch.fields
      for field in the_fields:
        name = field.name
        if name not in fields_dict.keys():
          fields_dict[name] = []


      # iterate over the mesh
      for cid in ch.mesh.cells.keys():
        upperZ = ch.mesh.cells[cid].upperZ
        lowerZ = ch.mesh.cells[cid].lowerZ
        upperX = (dx/2.0 + x_offset)
        lowerX = (-dx/2.0 + x_offset)
        upperY = (dy/2.0 + y_offset)
        lowerY = (-dy/2.0 + y_offset)

        origin = (x_offset, y_offset, (lowerZ + upperZ) / 2.0)

        local_pts = self._make_paraview_rect_prism((lowerX, upperX), (lowerY, upperY), (lowerZ, upperZ), origin)

        indices = list(range(point_index, point_index + 8))
        point_index += 8
        points.extend(local_pts)
        cells.extend(indices)
        offsets.append(len(cells))
        cell_types.append(12) # vtk hexahedron

        for field in the_fields:
          fields_dict[field.name].append(field.T[cid])
        if writeConductors:
          if len(ch.conductors) == ch.nZones:
            fields_dict['T_c_mean']+= [float(ch.conductors[cid].get_mean_temp())]

    # ---- Build XML tree ----
    vtk = ET.Element("VTKFile", type="UnstructuredGrid", version="0.1", byte_order="LittleEndian")
    ugrid = ET.SubElement(vtk, "UnstructuredGrid")
    piece = ET.SubElement(ugrid, "Piece",
                          NumberOfPoints=str(len(points)),
                          NumberOfCells=str(len(cell_types)))

    # Points
    pts_elem = ET.SubElement(piece, "Points")
    data = ET.SubElement(pts_elem, "DataArray",
                        type="Float64", NumberOfComponents="3", format="ascii")
    data.text = " ".join(f"{x} {y} {z}" for (x, y, z) in points)

    # Cells
    cells_elem = ET.SubElement(piece, "Cells")

    conn = ET.SubElement(cells_elem, "DataArray", type="Int32", Name="connectivity", format="ascii")
    conn.text = " ".join(str(i) for i in cells)

    offs = ET.SubElement(cells_elem, "DataArray", type="Int32", Name="offsets", format="ascii")
    offs.text = " ".join(str(o) for o in offsets)

    types = ET.SubElement(cells_elem, "DataArray", type="UInt8", Name="types", format="ascii")
    types.text = " ".join(str(t) for t in cell_types)

    # Cell data
    celldata = ET.SubElement(piece, "CellData", Scalars="field")

    for name in fields_dict.keys():
      farr = ET.SubElement(celldata, "DataArray",
                          type="Float64", Name=name, format="ascii")
      farr.text = " ".join(str(v) for v in fields_dict[name])

    # Write it
    tree = ET.ElementTree(vtk)
    ET.indent(tree, space="  ", level=0)
    tree.write(self.name, encoding="utf-8", xml_declaration=True)

  def _rotate_point_45(self, point: tuple[float, float, float], origin: tuple[float, float, float]) -> tuple[float, float, float]:
    """
    Rotates a point 45 degrees counter-clockwise around the z-axis.
    point: A tuple (x, y, z) representing the point to rotate.
    origin: A tuple (ox, oy, oz) representing the center of rotation.
    Returns a new tuple (x', y', z') with the rotated coordinates.
    """
    x, y, z = point
    ox, oy, oz = origin

    # Translate point to origin, rotate, then translate back
    sqrt2 = math.sqrt(2)
    translated_x = x - ox
    translated_y = y - oy
    new_x = (translated_x - translated_y) / sqrt2 + ox
    new_y = (translated_x + translated_y) / sqrt2 + oy
    return (new_x, new_y, z) # z does not change for rotation around z-axis

  def _make_paraview_rect_prism(self, x: tuple, y: tuple, z: tuple, origin: tuple[float, float, float]) -> str:
    """
    Makes a list of tuples of xyz points for
    a paraview rect prism.

    x : tuple
      -x and +x
    y : tuple
      -y and +y
    z : tuple
      -z and +z
    origin: tuple
      The (x, y, z) origin around which to rotate the prism.
    """
    lx, ux = x
    ly, uy = y
    lz, uz = z

    return [
      self._rotate_point_45((lx, ly, lz), origin),
      self._rotate_point_45((ux, ly, lz), origin),
      self._rotate_point_45((ux, uy, lz), origin),
      self._rotate_point_45((lx, uy, lz), origin),
      self._rotate_point_45((lx, ly, uz), origin),
      self._rotate_point_45((ux, ly, uz), origin),
      self._rotate_point_45((ux, uy, uz), origin),
      self._rotate_point_45((lx, uy, uz), origin),
    ]

    return the_rect_prism_str
