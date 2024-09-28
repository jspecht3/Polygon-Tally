from polygon_tally import *
import threading

import matplotlib
matplotlib.use("agg")

# global variables
path = "output/different-polygon-function/"

def write_output(string, file_path):
    f = open(file_path, 'a')
    f.write(string + '\n')
    f.close()

## parameters
order = 35
radius = 1
mesh_size = 5_000
side_num = 5

def function(x, y):
    return 10 * y**2 * np.cos(5 * x) - 10 * x * y

## creating the save path
loc_path = path + "{}sides/".format(side_num)
if not os.path.isdir(loc_path):
    os.mkdir(loc_path)
    os.mkdir(loc_path + "z_basis/")
    os.mkdir(loc_path + "k_basis/")


# plotting the analytical function
def plot_ana():
    k1 = ZApprox(side_num, radius, mesh_size)
    x, y = k1.x, k1.y
    
    # analytical
    Fa = function(x, y)
    bounds = [np.min(Fa), np.max(Fa)]
    k1.plotter(Fa, bounds, "Analytical", loc_path + "ana.png")


# calculating the values associated with the k basis
def calc_kbasis():
    print("K: starting {} sided polygon".format(side_num))
    k1 = ZApprox(side_num, radius, mesh_size)
    x, y = k1.x, k1.y
    
    # analytical
    Fa = function(x, y)
    bounds = [np.min(Fa), np.max(Fa)]
    k1.plotter(Fa, bounds, "Analytical", loc_path + "ana.png")
    
    # norms
    Fa_l2_norm = la.norm(Fa, 2)
    Fa_linf_norm = la.norm(Fa, np.inf)

    # approx
    Fk = 0
    
    
    for n in range(order + 1):
        print("K: --- n ---", n)
        for m in np.arange(-n, n + 1, 2):
            print("K: m", m)
            # k basis
            bk = KBasis(n, m, side_num, radius)
            bk._gen_all_cd(mesh_size)

            # ck
            ck = bk.num_ck_nm(mesh_size, function)
            ck_out = (str(n) + "," + str(m) + "," + str(ck))
            write_output(ck_out, loc_path + "cks.txt")
            
            Fk += np.float64(ck * bk.k_nm(x, y))

        print("K: norms")
        Dk = Fa - Fk
        l2_k = la.norm(Dk, 2) / Fa_l2_norm
        linf_k = la.norm(Dk, np.inf) / Fa_linf_norm

        print("K: text out")
        write_output(str(l2_k), loc_path + "l2k.txt")
        write_output(str(linf_k), loc_path + "linfk.txt")

        print("K: imgs out")
        k1.plotter(Fk, bounds, "K Order = {}".format(n), loc_path + "k_basis/k{}.png".format(n))


# calculating the values associated with the z basis
def calc_zbasis():
    print("Z: starting {} sided polygon".format(side_num))
    # setting approximations
    z1 = ZApprox(side_num, radius, mesh_size)
    x, y = z1.x, z1.y
    
    # analytical
    Fa = function(x, y)
    bounds = [np.min(Fa), np.max(Fa)]
        
    # norms
    Fa_l2_norm = la.norm(Fa, 2)
    Fa_linf_norm = la.norm(Fa, np.inf)

    # approx
    Fz = 0

    for n in range(order + 1):
        print("Z: --- n ---", n)
        for m in np.arange(-n, n + 1, 2):
            print("Z: m", m)
            # z basis
            bz = ZBasis(n, m, side_num, radius)
            bz._gen_all_cd(mesh_size)

            # cz
            cz = bz.num_cz_nm(mesh_size, function)
            cz_out = (str(n) + "," + str(m) + "," + str(cz))
            write_output(cz_out, loc_path + "czs.txt")
            
            Fz += np.float64(cz * bz.z_nm(x, y))

        print("Z: norms")
        Dz = Fa - Fz        
        l2_z = la.norm(Dz, 2) / Fa_l2_norm
        linf_z = la.norm(Dz, np.inf) / Fa_linf_norm
        
        print("Z: text out")
        write_output(str(l2_z), loc_path + "l2z.txt")
        write_output(str(linf_z), loc_path + "linfz.txt")

        print("Z: imgs out")
        z1.plotter(Fz, bounds, "Z Order = {}".format(n), loc_path + "z_basis/z{}.png".format(n))


# multithreading
if __name__ == "__main__":
    t1 = threading.Thread(target=calc_kbasis)
    t2 = threading.Thread(target=calc_zbasis)
    t3 = threading.Thread(target=plot_ana)

    t1.start()
    t2.start()
    t3.start()

    t1.join()
    t2.join()
    t3.join()

    print("Done!")
