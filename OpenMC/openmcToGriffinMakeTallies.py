import openmc

def make_xs_tallies(mesh, geometry, groups, mgxs_types, num_delayed_groups, cell_domains, uni_domains, legendre, fine_groups, tallies_file, doMesh, doCells, doUnis):
    #### PART 0 DECLARE TALLIES
    energy_groups = openmc.mgxs.EnergyGroups(groups)


    #### PART 1 DO CELL TALLIES
    if doCells:
      ## transport xs
      mgxs_cell_lib_transport = openmc.mgxs.Library(geometry)
      mgxs_cell_lib_transport.energy_groups = openmc.mgxs.EnergyGroups(fine_groups)
      mgxs_cell_lib_transport.mgxs_types = ["transport"]
      mgxs_cell_lib_transport.num_delayed_groups = num_delayed_groups
      mgxs_cell_lib_transport.domain_type = 'cell'
      mgxs_cell_lib_transport.domains = cell_domains
      mgxs_cell_lib_transport.legendre_order = 1
      mgxs_cell_lib_transport.by_nuclide = False
      mgxs_cell_lib_transport.build_library()
      mgxs_cell_lib_transport.add_to_tallies_file(tallies_file, merge=True)

      ## all other XS.
      mgxs_cell_lib = openmc.mgxs.Library(geometry)
      mgxs_cell_lib.energy_groups = energy_groups
      mgxs_cell_lib.mgxs_types = mgxs_types
      mgxs_cell_lib.num_delayed_groups = num_delayed_groups
      mgxs_cell_lib.domain_type = 'cell'
      mgxs_cell_lib.domains = cell_domains
      mgxs_cell_lib.legendre_order = legendre
      mgxs_cell_lib.by_nuclide = False
      # mgxs_cell_lib.check_library_for_openmc_mgxs()
      mgxs_cell_lib.build_library()
      mgxs_cell_lib.add_to_tallies_file(tallies_file, merge=True)
    else:
      mgxs_cell_lib_transport = []
      mgxs_cell_lib = []


    #### PART 2 DO MESH TALLIES
    if doMesh:
      ## Mesh tally for normal cross sections
      mgxs_mesh_lib = openmc.mgxs.Library(geometry)
      mgxs_mesh_lib.energy_groups = openmc.mgxs.EnergyGroups(groups)
      mgxs_mesh_lib.mgxs_types = mgxs_types
      mgxs_mesh_lib.num_delayed_groups = num_delayed_groups
      mgxs_mesh_lib.domain_type = 'mesh'
      mgxs_mesh_lib.domains = [mesh]
      mgxs_mesh_lib.legendre_order = legendre
      mgxs_mesh_lib.by_nuclide = False
      mgxs_mesh_lib.build_library()
      mgxs_mesh_lib.add_to_tallies_file(tallies_file, merge=True)

      ## need to get transport xs tallied in fine grid over our mesh
      mgxs_mesh_lib_transport = openmc.mgxs.Library(geometry)
      mgxs_mesh_lib_transport.energy_groups = openmc.mgxs.EnergyGroups(fine_groups)
      mgxs_mesh_lib_transport.mgxs_types = ["transport"]
      mgxs_mesh_lib_transport.num_delayed_groups = num_delayed_groups
      mgxs_mesh_lib_transport.domain_type = 'mesh'
      mgxs_mesh_lib_transport.domains = [mesh]
      mgxs_mesh_lib_transport.legendre_order = 1
      mgxs_mesh_lib_transport.by_nuclide = False
      mgxs_mesh_lib_transport.build_library()
      mgxs_mesh_lib_transport.add_to_tallies_file(tallies_file, merge=True)
    else:
      mgxs_mesh_lib_transport = []
      mgxs_mesh_lib = []



    #### PART 3 DO UNIVERSE TALLIES
    if doUnis:
      ## transport xs
      mgxs_uni_lib_transport = openmc.mgxs.Library(geometry)
      mgxs_uni_lib_transport.energy_groups = openmc.mgxs.EnergyGroups(fine_groups)
      mgxs_uni_lib_transport.mgxs_types = ["transport"]
      mgxs_uni_lib_transport.num_delayed_groups = num_delayed_groups
      mgxs_uni_lib_transport.domain_type = 'universe'
      mgxs_uni_lib_transport.domains = uni_domains
      mgxs_uni_lib_transport.legendre_order = 1
      mgxs_uni_lib_transport.by_nuclide = False
      mgxs_uni_lib_transport.build_library()
      mgxs_uni_lib_transport.add_to_tallies_file(tallies_file, merge=True)

      ## all other XS.
      mgxs_uni_lib = openmc.mgxs.Library(geometry)
      mgxs_uni_lib.energy_groups = energy_groups
      mgxs_uni_lib.mgxs_types = mgxs_types
      mgxs_uni_lib.num_delayed_groups = num_delayed_groups
      mgxs_uni_lib.domain_type = 'universe'
      mgxs_uni_lib.domains = uni_domains
      mgxs_uni_lib.legendre_order = legendre
      mgxs_uni_lib.by_nuclide = False
      # mgxs_uni_lib.check_library_for_openmc_mgxs()
      mgxs_uni_lib.build_library()
      mgxs_uni_lib.add_to_tallies_file(tallies_file, merge=True)

    else:
      mgxs_uni_lib_transport = []
      mgxs_uni_lib = []


    return tallies_file, mgxs_mesh_lib, mgxs_mesh_lib_transport, mgxs_cell_lib, mgxs_cell_lib_transport, mgxs_uni_lib, mgxs_uni_lib_transport
