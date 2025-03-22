from polygon_tally import *
import numpy as np
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt


def write_output(string, file_path):
    f = open(file_path, 'a')
    f.write(string + '\n')
    f.close()


path = "/home/joe/projects/polygon-tally/presentations/plots/"

# parameters
order = 35
radius = 1
mesh_size = 5000  # make 5k
side_num = 20

# making dirs
loc_path = path + f"{side_num}sides/"
if not os.path.isdir(loc_path):
    os.mkdir(loc_path)
    os.mkdir(loc_path + "z_basis/")
    os.mkdir(loc_path + "k_basis/")

# coordinates
z1 = ZApprox(side_num, radius, mesh_size)
k1 = KApprox(side_num, radius, mesh_size)
x, y = z1.x, z1.y

# analytical
Fa = base_input(x, y)

Fa_l2_norm = la.norm(Fa, 2)
Fa_linf_norm = la.norm(Fa, np.inf)


# bounds
cmap_range = 1.1 * (np.max(Fa) - np.min(Fa))
cmap_center = (np.max(Fa) + np.min(Fa)) / 2

lb = cmap_center - cmap_range / 2
ub = cmap_center + cmap_range / 2

# color map
cmap = mpl.colormaps['viridis']
norm = mpl.colors.Normalize(lb, ub)
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])


def plotter(z, title='', savepath=''):
    fig, ax = plt.subplots()
    plot = ax.contourf(x, y, z, norm=norm, levels=100)

    ax.axis('equal')
    cbar = plt.colorbar(
        mappable=sm,
        ticks=np.linspace(int(np.min(Fa)), int(np.max(Fa)), 7),
        ax=ax
    )

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.spines[['top', 'bottom', 'left', 'right']].set_visible(False)

    if title:
        plt.title=title
    if savepath:
        plt.savefig(f"{savepath}.png", dpi=600)
    plt.close()


# do the analytical save
plotter(Fa, title="Analytical", savepath=loc_path + "ana")

Fz = 0
Fk = 0

# coefficient calcs
for n in range (order + 1):
    print(f"n = {n}")
    for m in np.arange(-n, n + 1, 2):
        #print(f"m = {m}")
        # z basis
        bz = ZBasis(n, m, side_num, radius)
        bz._gen_all_cd(mesh_size)

        cz = bz.num_cz_nm(mesh_size)
        cz_out = (str(n) + "," + str(m) + "," + str(cz))
        write_output(cz_out, loc_path + "czs.txt")

        Fz += np.float64(cz * bz.z_nm(x, y))

        # k basis
        bk = KBasis(n, m, side_num, radius)
        bk._gen_all_cd(mesh_size)

        ck = bk.num_ck_nm(mesh_size)
        ck_out = (str(n) + "," + str(m) + "," + str(ck))
        write_output(ck_out, loc_path + "cks.txt")

        Fk += np.float64(ck * bk.k_nm(x, y))

    print("norms")
    # norm calcs
    Dz = Fa - Fz
    Dk = Fa - Fk

    l2_z = la.norm(Dz, 2) / Fa_l2_norm
    l2_k = la.norm(Dk, 2) / Fa_l2_norm

    linf_z = la.norm(Dz, np.inf) / Fa_linf_norm
    linf_k = la.norm(Dk, np.inf) / Fa_linf_norm

    print("text out")
    # writing text outputs
    write_output(str(l2_z), loc_path + "l2z.txt")
    write_output(str(l2_k), loc_path + "l2k.txt")
    write_output(str(linf_z), loc_path + "linfz.txt")
    write_output(str(linf_k), loc_path + "linfk.txt")

    print("imgs out")
    # img outputs
    plotter(Fz, title=f"Z Order = {n}", savepath=loc_path + f"z_basis/z{n}")
    plotter(Fk, title=f"K Order = {n}", savepath=loc_path + f"k_basis/k{n}")
