# packages
import numpy as np
from numpy import linalg as la
from numpy import arctan as atan
from numpy import sin, cos, tan, pi

import matplotlib as mpl
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.ticker import LinearLocator

import math
from math import factorial as fact

import csv
import time
import scipy as sy
import pandas as pd
from PIL import Image

import warnings
warnings.filterwarnings('ignore')

np.set_printoptions(edgeitems=30, linewidth=100000, 
    formatter=dict(float=lambda x: "%.5g" % x))

from scipy.integrate import quad, dblquad, nquad
import random
import pytest

# zernikes
def m_checker(n,m):
    ms = np.arange(-n,n+1,2)
    if (m in ms) == False:
        raise ValueError('invalid m')
    return

def n_mn(n,m):
    if m != 0:
        return (2*(n+1))**(1/2)
    else:
        return (n+1)**(1/2)

def r_mn(n,m,rho):
    stop = int((n-abs(m))/2)
    toRet = 0
    for k in range(stop+1):
        top = (-1)**k * fact(int(n-k)) * rho**(n-2*k)
        bot = fact(k) * fact(int((n+abs(m))/2 - k)) * fact(int((n-abs(m))/2 - k))
        toRet += top/bot
    return toRet

def zernike_m(n,m,rho,phi):
    m_checker(n,m)
    if m >= 0:
        return n_mn(n,m) * r_mn(n,m,rho) * cos(m*phi)
    else:
        return -1 * n_mn(n,m) * r_mn(n,m,rho) * sin(m*phi)


# coordinate generators
def gen_rhophi(size):
    rho = np.linspace(0,1,size)
    phi = np.linspace(0,2*pi,size)
    rho,phi = np.meshgrid(rho,phi)
    return rho,phi

def gen_uv(rho,phi):
    u = rho * cos(phi)
    v = rho * sin(phi)
    return u,v

def r_alpha(phi,p,R0):
    alpha = pi / p
    
    def _u_alpha(phi):
        x = (phi + alpha) / (2*alpha)
        x = x.astype(int)
        return phi - x*(2*alpha)

    U = _u_alpha(phi)
    return (R0 * cos(alpha)) / cos(U)

def gen_rtheta(rho,phi,p,R0):
    R = r_alpha(phi,p,R0)
    r = rho * R
    theta = phi
    return r,theta

def gen_xy(rho,phi,p,R0):
    r,theta = gen_rtheta(rho,phi,p,R0)
    x = r * cos(theta)
    y = r * sin(theta)
    return x,y

def gen_polar_differentials(r,theta):
    dr = r[:,1]
    _,dr = np.meshgrid(dr,dr)
    dtheta = theta[1,1]
    return dr,dtheta


# plotting
def plotter(x,y,z,title='',save_fig=False,save_path='',cbar_lim=True):
    fig,ax = plt.subplots()
    ax.set_title(title)
    
    plot = ax.contourf(x,y,z)
    cbar = fig.colorbar(plot)
    if cbar_lim == True: plot.set_clim(-300,180)

    if save_fig == True: plt.savefig(save_path)
    plt.show()

def plot_twinx(data1,data2,title,xlabel,ylabel1,ylabel2,go_log_mode=False,save_fig=False,save_path=''):
    fig,ax1 = plt.subplots()

    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    
    ax1.set_ylabel(ylabel1, color='r')
    ax1.plot(data1, color='r')
    ax1.tick_params(axis='y', labelcolor='r')
    
    ax2 = ax1.twinx()
    ax2.set_ylabel(ylabel2, color='b')
    ax2.plot(data2, color='b')
    ax2.tick_params(axis='y', labelcolor='b')

    if go_log_mode == True:
        ax1.set_yscale('log')
        ax2.set_yscale('log')

    if save_fig == True:
        plt.savefig(save_path, dpi=600)

    plt.show()

def polygon_area(p,R0): return p * R0**2 * cos(pi/p)**2 * tan(pi/p)


# transformation functions
def base_input(x,y): return 2*x**2 - y**2 + x**2*y - 4*x*y**2 + 5*x*y - 3*x + 5*y

## K
# general function
def k_m(n,m,r,theta,p,R0):
    phi = theta
    R = r_alpha(phi,p,R0)
    rho = r / R

    K = zernike_m(n,m,rho,phi)
    return K

# weighting factor
def ck(n,m,r,theta,p,R0,f=base_input):
    x = r * cos(theta)
    y = r * sin(theta)

    F = f(x,y)
    K = k_m(n,m,r,theta,p,R0)

    #Integration Parameters
    R = r_alpha(theta,p,R0)
    
    dr = r[:,1]
    _,dr = np.meshgrid(dr,dr)
    dtheta = theta[1,1]

    dmu = r * dr * dtheta / R**2 / pi

    return np.sum(K * F * dmu)

# basis vector
def fk_m(n,m,r,theta,p,R0,f=base_input):
    c = ck(n,m,r,theta,p,R0,f)
    K = k_m(n,m,r,theta,p,R0)
    return c*K

def fk_n(n,r,theta,p,R0,f=base_input):
    ms = np.arange(-n,n+1,2)
    F = 0
    for m in ms:
        F += fk_m(n,m,r,theta,p,R0,f)
    return F

def fk(n,r,theta,p,R0,f=base_input):
    F = 0
    for i in range(n+1):
        F += fk_n(i,r,theta,p,R0,f)
    return F


## H
# general function
def h_m(n,m,r,theta,p,R0):
    phi = theta
    R = r_alpha(phi,p,R0)
    rho = r / R

    H = zernike_m(n,m,rho,phi) / R
    return H

# weighting factor
def ch(n,m,r,theta,p,R0,f=base_input):
    x = r * cos(theta)
    y = r * sin(theta)

    F = f(x,y)
    H = h_m(n,m,r,theta,p,R0)

    #Integration Parameters
    R = r_alpha(theta,p,R0)
    
    dr = r[:,1]
    _,dr = np.meshgrid(dr,dr)
    dtheta = theta[1,1]

    dmu = r * dr * dtheta / pi

    return np.sum(H * F * dmu)

# basis vector
def fh_m(n,m,r,theta,p,R0,f=base_input):
    c = ch(n,m,r,theta,p,R0,f)
    H = h_m(n,m,r,theta,p,R0)
    return c*H

def fh_n(n,r,theta,p,R0,f=base_input):
    ms = np.arange(-n,n+1,2)
    F = 0
    for m in ms:
        F += fh_m(n,m,r,theta,p,R0,f)
    return F

def fh(n,r,theta,p,R0,f=base_input):
    F = 0
    for i in range(n+1):
        F += fh_n(i,r,theta,p,R0,f)
    return F
