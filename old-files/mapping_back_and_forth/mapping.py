from glip_glorp import *
SAVE_FIGS = True

# initial polar definition
rho = np.linspace(0,1,1000)
phi = np.linspace(0,2*pi,1000)
rho,phi = np.meshgrid(rho,phi)

# polar -> cartesian
u = rho * cos(phi)
v = rho * sin(phi)

# cartesian -> polar
r = (u**2 + v**2)**(1/2)
theta = np.arctan2(v,u)

# polar -> cartesian, 2
x = r * cos(theta)
y = r * sin(theta)

# plotting
fig,ax = plt.subplots()
plot = ax.contourf(u,v,0*u)
cbar = fig.colorbar(plot)
plt.title("Initial Cartesian")
if SAVE_FIGS: plt.savefig("initial", dpi=600)

fig,ax = plt.subplots()
plot = ax.contourf(x,y,0*x)
cbar = fig.colorbar(plot)
plt.title("Transformed Cartesian")
if SAVE_FIGS: plt.savefig("transformed", dpi=600)
