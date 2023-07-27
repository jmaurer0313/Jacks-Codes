#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 11:38:24 2018

Script to calculate the interaction energy of two extended dipoles

@author: jkittell
"""

# Cartesian vector system: unrotated dipoles anti-aligned
#    x
#    |
#  -\+/- y
#    |
#
# Extended dipoles: vectors with a plus charge at the head and a minus charge
#   at the tail.

import math as ma
Pi = ma.pi
import numpy as np
from numpy import linalg as la



def mag(x): 
    return ma.sqrt(sum(i**2 for i in x))

# Debye:
D = 3.33e-30
mumonomer = 12.8 * D

# charge
c = .38*1.602176e-19

# separation distance: (aligned with z)
#sep = 6e-10
sep = 7.7e-10

#vector len:

l = 7e-10
#l = 7e-10


# dipole orientations:
def dip_vecs(phi, theta):

    return [np.array([ma.cos(phi/2)*ma.cos(theta/2), 
                      -ma.sin(phi/2), 
                      -ma.cos(phi/2)*ma.sin(theta/2)]),
            np.array([ma.cos(-phi/2)*ma.cos(-theta/2), 
                      -ma.sin(-phi/2), 
                      -ma.cos(-phi/2)*ma.sin(-theta/2)])]

    
def J_exten(phi_r, theta_r, sep, l, q):
    '''Coupling for extended dipoles:
    rads, rads, meters, meters, coulombs -> Joules'''
    k = 8.99e9
    
    #determine vector locations of 4 charges:
    vecs = dip_vecs(phi_r, theta_r)
    head1 = np.array([0,0,sep/2]) + vecs[0]*(l/2)
    head2 = np.array([0,0,-sep/2]) + vecs[1]*(l/2)
    tail1 = np.array([0,0,sep/2]) - vecs[0]*(l/2)
    tail2 = np.array([0,0,-sep/2]) - vecs[1]*(l/2)
    
    #calculate energy:
    head_head = k * q**2 / mag(head1-head2)
    tail_tail = k * q**2 / mag(tail1-tail2)
    head_tail =-k * q**2 / mag(head1-tail2)
    
    #why not the factor of J2NuBar Here?
    return (head_head + tail_tail + 2 * head_tail)

#
#test = []
#thetas = range(0, 105, 15)
#phis = range(0, 105, 15)
#Js = []
#
#for theta in thetas:
#    for phi in phis:
#        Js.append(J_exten(phi,theta,sep, l, c))


def J_point(phi_r, theta_r, sep, l, q):
    "coupling for point dipoles"
    k = 8.99e9
    mumonomer = l * q

    R12 = np.array([0,0,1]) * sep
    mu1, mu2 = dip_vecs(phi_r, theta_r)
    
    mu1 = mu1 * mumonomer
    mu2 = mu2 * mumonomer
    
    return (k/(la.norm(R12)**3))*(np.dot(mu1, mu2)- \
                                            3*(np.dot(R12,mu1)*(np.dot(mu2,R12)/(la.norm(R12))**2))) 
#TEST# print(J_point(0,0,sep))


####  Plotting below....  #####

#'''
#### Phi Plots ###
#phi_exten = []
#for i in range(360):
#    phi_exten.append([i, J_exten(i*(Pi/180),0,sep,l,c)])
#phi_exten = np.array(phi_exten)
#
#
#phi_point = []
#for i in range(360):
#    phi_point.append([i, J_point(i*(Pi/180), 0, sep, l, c)])
#phi_point = np.array(phi_point)
#
#### Theta Plots ###
#
#theta_exten = []
#for i in range(360):
#    theta_exten.append([i, J_exten(0,i*(Pi/180),sep,l,c)])
#theta_exten = np.array(theta_exten)
#
#
#theta_point = []
#for i in range(360):
#    theta_point.append([i, J_point(0, i*(Pi/180), sep, l, c)])
#theta_point = np.array(theta_point)

##  3D Plots  #
#
#dat_exten = np.zeros((180,360))
#dat_point = np.zeros((180,360))
#
#for i in range(180):
#    for j in range(360):
#        dat_exten[i, j] = J_exten(i, j, sep, l, c)
#        dat_point[i, j] = J_point(i, j, sep, l, c)
#

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as A3D
from matplotlib import cm

'''

# 2D Plots #

### Tilt Only and Twist Only ####
plt.figure(1, figsize = (7,5))
print(len(phi_exten))
plt.plot(phi_exten[:,0], phi_exten[:,1])
plt.plot(phi_point[:,0], phi_point[:,1])
plt.xlabel("Phi (Degrees)")
plt.ylabel("J (Joules)")
plt.legend(["Extended", "Point"])
plt.title("Coupling Strength: Twist Only")

plt.axis((0,360,-7e-20,7e-20))
plt.figure(2, figsize = (7,5))
plt.plot(theta_exten[:,0], theta_exten[:,1])
plt.plot(theta_point[:,0], theta_point[:,1])
plt.xlabel("Theta (Degrees)")
plt.ylabel("J (Joules)")
plt.legend(["Extended", "Point"])
plt.title("Coupling Strength: Tilt Only")
plt.axis((0,360,0,40e-20))
plt.show()
### Gif Creation ###

#alphalist = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
#file_names = []
#for i in range(180):
#    file_names.append(str(i+1)+"_")
#angle = 0
#for file_name in file_names:
#    



#### 3D Plots ####


#plt.figure(3, figsize = (10,10))
#plt.contour(dat_exten, levels = 15)
#plt.figure(4, figsize = (10,8))
#plt.contour(dat_point, levels = 15)


#
#fig = plt.figure(5, figsize = (13, 10))
#ax1 = fig.add_subplot(111, projection='3d')
#ax1.set_xlabel('Phi (Degrees)')
#ax1.set_ylabel('Theta (Degrees)')
#ax1.set_zlabel('J (Joules)')
#X, Y = np.meshgrid(np.arange(360), np.arange(180))
#Z = np.array(dat_exten)
#ax1.plot_wireframe(X, Y, Z, ccount = 0)
#Z = np.array(dat_point)
#ax1.plot_wireframe(X, Y, Z, ccount = 0, color = "G")
#ax1.set_zlim(-.5e-19, 2e-19)
#for half_ang in range(180):
#    ax1.view_init(8, 2*half_ang)
#    fig.savefig("/Users/jkittell/Dropbox/Justin Project/Modeling/Extended Dipole/" + str(half_ang) + "_b.png")


#ax2 = fig.add_subplot(122, projection='3d')
#X, Y = np.meshgrid(np.arange(360), np.arange(180))
#Z = np.array(dat_point)
#ax2.plot_wireframe(X, Y, Z, ccount = 0)
#ax2.set_zlim(-.5e-19, 2e-19)
#ax2.view_init(15, 90)
#
#plt.subplots_adjust(wspace = 0.001)
#

#fig2 = plt.figure(6, figsize = (10, 4))
#ax1 = fig2.add_subplot(131)
#X, Y = np.meshgrid(np.arange(360), np.arange(180))
#Z = np.array(dat_exten)
##ax1.contour(X, Y, Z, levels = np.arange(-100, 400)*1e-19 / 100)
#ax1.pcolor(Z, vmin = -2e-19, vmax = 2e-19, cmap = "seismic")
#
#ax2 = fig2.add_subplot(132)
#X, Y = np.meshgrid(np.arange(360), np.arange(180))
#Z = np.array(dat_point)
#ax2.pcolor(Z, vmin = -2e-19, vmax = 2e-19, cmap = "seismic")
#
#ax3 = fig2.add_subplot(133)
#X, Y = np.meshgrid(np.arange(360), np.arange(180))
#Z = np.array(dat_exten)-np.array(dat_point)
#ax3.pcolor(Z, vmin = -2e-19, vmax = 2e-19, cmap = "seismic")
#
#plt.tight_layout()
#plt.show()
'''