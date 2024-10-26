import openmc
import matplotlib.pyplot as plt
import numpy as np
import csv

# materials
## core, weapons grade Pu
core = openmc.Material()

enrichment_tru = 15.5 / 100
enrichment_u = 1 - enrichment_tru

core.add_nuclide("U234", 0.001 * enrichment_u)
core.add_nuclide("U235", 0.2 * enrichment_u)
core.add_nuclide("U238", 99.8 * enrichment_u)

core.add_nuclide("Pu238", 0.01 * enrichment_tru)
core.add_nuclide("Pu239", 93.81 * enrichment_tru)
core.add_nuclide("Pu240", 5.81 * enrichment_tru)
core.add_nuclide("Pu241", 0.35 * enrichment_tru)
core.add_nuclide("Pu242", 0.02 * enrichment_tru)

core.set_density("g/cc", 20)

## fill
na = openmc.Material()
na.add_element("Na", 100)
na.set_density("g/cc", 0.971)

## HT9, clad
ht9 = openmc.Material()

ht9.add_element("C", 0.16)
ht9.add_element("Si", 0.04)
ht9.add_element("Mn", 0.58)
ht9.add_element("Cr", 12.20)
ht9.add_element("Mo", 0.90)
ht9.add_element("W", 0.50)
ht9.add_element("V", 0.29)
ht9.add_element("Ni", 0.69)
ht9.add_element("S", 0.002)
ht9.add_element("P", 0.003)
ht9.add_element("N", 0.106)
ht9.add_element("Fe", 84.529)

ht9.set_density("g/cc", 6.551)

## helium
he = openmc.Material()
he.add_element("He", 1)

he.set_density("g/cc", 0.166 / 1000)

## xml
materials = openmc.Materials([core, na, ht9, he])
#materials.export_to_xml()

## color by
colors = {
        core : 'salmon',
        na : 'steelblue',
        ht9 : 'darkgray',
        he : 'red',
}


# geometry
## fuel pin
### parameters
pellet_diameter = 5.5685 / 10

clad_di = 6.43 / 10
clad_do = 7.55 / 10

pin_cell_f2f = 8.9074 / 10
pin_pitch = 1.180 * pellet_diameter
duct_f2f = 15.710

### cells
fuel_or = openmc.ZCylinder(r = pellet_diameter / 2)
clad_ir = openmc.ZCylinder(r = clad_di / 2)
clad_or = openmc.ZCylinder(r = clad_do / 2)
pin_hex = openmc.model.HexagonalPrism(
        edge_length = pin_cell_f2f * 3**(-0.5))

fuel_cell = openmc.Cell(fill = core, region = -fuel_or) 
gap_cell = openmc.Cell(fill = na, region = +fuel_or & -clad_ir)
clad_cell = openmc.Cell(fill = ht9, region = +clad_ir & -clad_or)
outer_cell = openmc.Cell(fill = na, region = +clad_or)

pin_cell = openmc.Cell(
        fill = openmc.Universe(
            cells = [fuel_cell,
                     gap_cell,
                     clad_cell,
                     outer_cell]))
pin_univ = openmc.Universe(cells = [pin_cell])

pin_univ.plot(color_by = "material", width = [1.1, 1.1], colors = colors)
plt.savefig("pin-univ.png", dpi=600)
plt.close()

## assembly
assembly_hex = openmc.HexLattice()

assembly_hex.center = (0,0)
assembly_hex.pitch = (pin_pitch * 3**(0.5),)
assembly_hex.outer = openmc.Universe(cells = [openmc.Cell(fill=he)])

layers = []

for i in reversed(range(1,10)):
    layers.append([pin_univ] * (6*i))
layers.append([pin_univ])

assembly_hex.universes = layers

### cell
polygon_radius = pin_pitch * 10 * 3**(0.5)
outer_hex = openmc.model.HexagonalPrism(
        edge_length = polygon_radius,
        orientation = 'y',
        boundary_type = "reflective")

assembly_cell = openmc.Cell(fill = assembly_hex, region = -outer_hex)

### plotting
assembly_univ = openmc.Universe(cells=[assembly_cell])
assembly_univ.plot(
        color_by = "material",
        width = [25,25],
        colors = colors,
        pixels = [500,500])
plt.savefig("assembly-univ.png", dpi=600)
plt.close()

geometry = openmc.Geometry([assembly_cell])
#geometry.export_to_xml()


# settings
settings = openmc.Settings()

source = openmc.stats.Box(
        (-polygon_radius,-polygon_radius,-polygon_radius), 
        (polygon_radius,polygon_radius,polygon_radius))
src = openmc.IndependentSource(space = source)
settings.soruce = [src]

settings.particles = 100
settings.batches = 500
settings.inactive = 250

#settings.export_to_xml()


# tallies
tallies = openmc.Tallies()

## mesh
mesh = openmc.RegularMesh()

mesh_len = 250
mesh_dimensions = (mesh_len, mesh_len)
mesh.dimension = mesh_dimensions

mesh_scope = polygon_radius
mesh.lower_left = (-mesh_scope, -mesh_scope)
mesh.upper_right = (mesh_scope, mesh_scope)

mesh_filter = openmc.MeshFilter(mesh)

## tally
mesh_tally = openmc.Tally()
mesh_tally.scores = ['kappa-fission', 'flux']
mesh_tally.filters = [mesh_filter]

tallies.append(mesh_tally)

## zernike
zernike_filter = openmc.ZernikeFilter(
        order=15,
        r = polygon_radius
)

zernike_tally = openmc.Tally()
zernike_tally.filters = [zernike_filter]
zernike_tally.scores = ['kappa-fission', 'flux']

tallies.append(zernike_tally)

# model
model = openmc.Model(
        geometry = geometry,
        materials = materials,
        settings = settings,
        tallies = tallies
)
model.export_to_xml()

# data needed in other files
with open('transfer.csv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter=' ')
    
    writer.writerow(["mesh_dimensions",
                     mesh_dimensions[0],
                     mesh_dimensions[1]])
    writer.writerow(["pin_f2f", pin_cell_f2f])
    writer.writerow(["polygon_radius", polygon_radius])
