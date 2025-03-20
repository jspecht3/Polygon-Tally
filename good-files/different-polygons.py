from polygon_tally import *
import os

path = "output/different-polygons/"

def write_output(string, file_path):
    f = open(file_path, 'a')
    f.write(string + '\n')
    f.close()

# parameters
order = 35
radius = 1
mesh_size = 5_000

for side_num in [3, 4, 6, 8, 20]:
    loc_path = path + "{}sides/".format(side_num)
    if not os.path.isdir(loc_path):
        os.mkdir(loc_path)
        os.mkdir(loc_path + "z_basis/")
        os.mkdir(loc_path + "k_basis/")

    z1 = ZApprox(side_num, radius, mesh_size)
    k1 = KApprox(side_num, radius, mesh_size)
    x, y = z1.x, z1.y

    print("starting {} sided polygon".format(side_num))
    # analytical
    Fa = base_input(x, y)
    bounds = [np.min(Fa), np.max(Fa)]

    Fa_l2_norm = la.norm(Fa, 2)
    Fa_linf_norm = la.norm(Fa, np.inf)

    # approx
    Fz = 0
    Fk = 0

    # plot ana
    k1.plotter(Fa, bounds, "Analytical", loc_path + "ana.png")

    # coefficient calcs
    for n in range (order + 1):
        print("n", n)
        for m in np.arange(-n, n + 1, 2):
            print("m", m)
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
        z1.plotter(Fz, bounds, "Z Order = {}".format(n), loc_path + "z_basis/z{}.png".format(n))
        k1.plotter(Fk, bounds, "K Order = {}".format(n), loc_path + "k_basis/k{}.png".format(n))
