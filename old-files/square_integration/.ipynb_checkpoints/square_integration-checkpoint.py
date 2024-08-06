from glip_glorp import *

# Global Variables
SAVE_FIGS = 1
p = 4
R0 = 1
alpha = pi / p
side_length = R0 / 2**(1/2)


### Begin U_alpha and R_alpha

# functions
u_alpha = lambda phi : phi - ((phi + alpha)/(2*alpha)).astype(int) * 2 * alpha
new_r_alpha = lambda phi : R0 * cos(alpha) / cos(u_alpha(phi))

## Plotting U_alpha and R_alpha
# variables for plotting
phi = np.linspace(0,2*pi,1000)
y_u = u_alpha(phi)
y_r = new_r_alpha(phi)

# plotting
plt.plot(phi, y_u, label='U')
plt.plot(phi, y_r, label='R')
[plt.axvline(drop, color='k', ls=(0,(5,3)), lw=1) for drop in (pi*np.array([1/4,3/4,5/4,7/4]))]
plt.legend()
#plt.title("Plot of U_alpha and R_alpha")
if SAVE_FIGS: plt.savefig('images/ur_alpha.png', dpi=600)


### End U_alpha and R_alpha

### Begin Zernike Functions

n_nm = lambda n,m : ((2*(n+1))/(1+(m==0)))**(1/2)

def r_nm(rho,n,m):
    stop = int((n-abs(m))/2)
    R_nm = 0

    for k in range(stop+1):
        top = (-1)**k * fact(n-k) * rho**(n-2*k)
        bot = fact(k) * fact(int((n+abs(m))/2 - k)) * fact(int((n-abs(m))/2 - k))
        R_nm += top/bot

    return R_nm

def ana_zernike(rho,phi,n,m):
    m_checker(n,m)
    if m >= 0:
        return n_nm(n,m) * r_nm(rho,n,m) * cos(m*phi)
    else:
        return -1* n_nm(n,m) * r_nm(rho,n,m) * sin(m*phi)

### End Zernike Functions

### Begin Basis Vectors

def k_nm(x,y,n,m):
    r = (x**2 + y**2)**(1/2)
    theta = np.arctan2(y,x) + (2*pi * (y < 0))
    
    R_alpha = new_r_alpha(theta)
    rho = r / R_alpha
 
    return ana_zernike(rho,theta,n,m)

def ck_nm(n,m,f=base_input):
    def _toIntegrate(x,y,n,m):
        K = k_nm(x,y,n,m)
        F = f(x,y)

        theta = np.arctan2(y,x) + (2*pi * (y < 0))
        R = new_r_alpha(theta)

        dmu = 1 / pi / R**2

        return K * F * dmu

    ck = nquad(_toIntegrate, [[-side_length,side_length],[-side_length,side_length]],args=(n,m))[0]
    
    return ck
