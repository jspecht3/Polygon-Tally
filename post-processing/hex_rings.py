import numpy as np
import matplotlib.pyplot as plt

pellet_diameter = 5.5685 / 10
pin_pitch = 1.180 * pellet_diameter
pitch = pin_pitch * 3**(0.5)
radius = pitch / 2 / np.cos(np.radians(30))


def next1(ring):
    x0, y0 = cs[ring][-1]       # previous center
    x1 = x0 + pitch * 3**(0.5) / 2  # new x
    y1 = y0 - pitch / 2         # new y
    cs[ring].append(np.array([x1, y1]))


def next2(ring):
    x0, y0 = cs[ring][-1]
    x1 = x0
    y1 = y0 - pitch
    cs[ring].append(np.array([x1, y1]))


def next3(ring):
    x0, y0 = cs[ring][-1]
    x1 = x0 - pitch * 3**(0.5) / 2
    y1 = y0 - pitch / 2
    cs[ring].append(np.array([x1, y1]))


def next4(ring):
    x0, y0 = cs[ring][-1]
    x1 = x0 - pitch * 3**(0.5) / 2
    y1 = y0 + pitch / 2
    cs[ring].append(np.array([x1, y1]))


def next5(ring):
    x0, y0 = cs[ring][-1]
    x1 = x0
    y1 = y0 + pitch
    cs[ring].append(np.array([x1, y1]))


def next6(ring):
    x0, y0 = cs[ring][-1]
    x1 = x0 + pitch * 3**(0.5) / 2
    y1 = y0 + pitch / 2
    cs[ring].append(np.array([x1, y1]))


def plotter():
    plt.figure(figsize=(6, 6))
    for ring in cs:
        for cell in cs[ring]:
            plt.scatter(cell[0], cell[1], color='b')
    plt.xlim(-10, 10)
    plt.ylim(-12, 12)
    plt.savefig('hex-lattice-proof.png', dpi=600)
    plt.show()


cs = {0: [np.array([0, 0])]}
for ring in range(1, 10):
    cs[ring] = [np.array([0, pitch * ring])]

    for cell in range(1, 6*ring):
        if cell <= ring * 1:
            next1(ring)
            continue
        if cell <= ring * 2:
            next2(ring)
            continue
        if cell <= ring * 3:
            next3(ring)
            continue
        if cell <= ring * 4:
            next4(ring)
            continue
        if cell <= ring * 5:
            next5(ring)
            continue
        if cell <= ring * 6:
            next6(ring)
            continue


def writer():
    with open("post-processor.txt", 'w') as f:
        f.write("[Postprocessors]\n")
        for ring in range(10):
            for cell_num in range(len(cs[ring])):
                cell_id = f"cell_{ring}_{cell_num}"
                cell_x, cell_y = cs[ring][cell_num]

                f.write(f"  [{cell_id}]\n" +
                        "    type = PointValue\n" +
                        "    variable = kappa_fission\n" +
                        f"    point = \'{cell_x} {cell_y} 0.0\'\n" +
                        "  []\n")
        f.write("[]")


# plotter()
# writer()
