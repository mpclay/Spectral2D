# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Copyright (C) 2014 by Matthew Clay <mpclay@gmail.com>
#
# File: TaylorGreen.py
# Author: Matthew Clay
#
# This script post-process Taylor-Green restart files to assess the temporal
# accuracy of the Spectral2D code. The Taylor-Green vortex is an exact solution
# of the Navier-Stokes equations on a periodic domain, and has solution:
#
#   psi(x,y,t) = sin(x)*sin(y)*F(t)
#   w(x,y,t) = 2*sin(x)*sin(y)*F(t)
#   u(x,y,t) = sin(x)*cos(y)*F(t)
#   v(x,y,t) = -cos(x)*sin(y)*F(t),
#
# where psi is the streamfunction, w is the vorticity, u is the velocity in the
# x direction, v is the velocity in the y direction, F(t) = exp(-2*nu*t), and
# nu is the kinematic viscosity.
#
# NOTE: In this analysis we assume that the final time and viscosity for all of
# the simulations are the same. Doing temporal analysis without these parameters
# being the same would be useless.
#

import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

# Function to calculate the L2 error between the exact and numerical solution.
def L2Error(dx, dy, psiE, wE, uE, vE, psiN, wN, uN, vN):
    psiDiff = np.subtract(psiE, psiN)
    wDiff = np.subtract(wE, wN)
    uDiff = np.subtract(uE, uN)
    vDiff = np.subtract(vE, vN)
    errPsi = np.sqrt(np.sum(np.square(psiDiff))*dx*dy)
    errW = np.sqrt(np.sum(np.square(wDiff))*dx*dy)
    errU = np.sqrt(np.sum(np.square(uDiff))*dx*dy)
    errV = np.sqrt(np.sum(np.square(vDiff))*dx*dy)
    return errPsi, errW, errU, errV

# Function to calculate the LInf error between the exact and numerical solution.
def LInfError(psiE, wE, uE, vE, psiN, wN, uN, vN):
    errPsi = np.max(np.abs(np.subtract(psiE, psiN)))
    errW = np.max(np.abs(np.subtract(wE, wN)))
    errU = np.max(np.abs(np.subtract(uE, uN)))
    errV = np.max(np.abs(np.subtract(vE, vN)))
    return errPsi, errW, errU, errV

# List of restart files to use when assessing temporal accuracy.
rests = ['REST_000001-00100.h5', \
         'REST_000001-01000.h5', \
         'REST_000001-10000.h5']
nRest = len(rests)

# Grid size and spacing.
fid = h5.File('GRID.h5', 'r')
nx = fid['Header'].attrs['nx']
ny = fid['Header'].attrs['ny']
dx = 2.0*np.pi/float(nx)
dy = 2.0*np.pi/float(ny)

# Get the grid.
x = np.zeros((nx, ny), np.float64)
y = np.zeros((nx, ny), np.float64)
x[:,:] = np.transpose(fid['Domain_00001']['x'])
y[:,:] = np.transpose(fid['Domain_00001']['y'])
fid.close()

# Form the necessary components to calculate the exact solution.
sinx = np.sin(x)
cosx = np.cos(x)
siny = np.sin(y)
cosy = np.cos(y)

# Calculate the exact solutions.
fid = h5.File(rests[0], 'r')
tEnd = fid['Header'].attrs['time']
nu = fid['Header'].attrs['nu']
F = np.exp(-2.0*nu*tEnd)
psiE = np.multiply(sinx, siny)*F
wE = 2.0*np.multiply(sinx, siny)*F
uE = np.multiply(sinx, cosy)*F
vE = -1.0*np.multiply(cosx, siny)*F
fid.close()

# Figure out the time steps and errors for each simulation.
dt = np.zeros(nRest, dtype=np.float64)
L2Err = np.zeros((nRest, 4), dtype=np.float64)
LInfErr = np.zeros((nRest, 4), dtype=np.float64)
psiN = np.zeros((nx, ny), np.float64)
wN = np.zeros((nx, ny), np.float64)
uN = np.zeros((nx, ny), np.float64)
vN = np.zeros((nx, ny), np.float64)
for i in range(nRest):
    # Get the time step for this solution.
    fid = h5.File(rests[i], 'r')
    tEnd = fid['Header'].attrs['time']
    nadv = fid['Header'].attrs['nadv']
    dt[i] = tEnd/float(fid['Header'].attrs['nadv'])
    #
    # Get the data for this simulation.
    psiN[:,:] = np.transpose(fid['FlowData']['Psi'])
    wN[:,:] = np.transpose(fid['FlowData']['Omega'])
    uN[:,:] = np.transpose(fid['FlowData']['u'])
    vN[:,:] = np.transpose(fid['FlowData']['v'])
    #
    # Calculate the L2 error for this simulation.
    psiL2, wL2, uL2, vL2 = L2Error(dx, dy, psiE, wE, uE, vE, psiN, wN, uN, vN)
    L2Err[i,0] = psiL2
    L2Err[i,1] = wL2
    L2Err[i,2] = uL2
    L2Err[i,3] = vL2
    #
    # Calculate the LInf error for the simulation.
    psiInf, wInf, uInf, vInf = LInfError(psiE, wE, uE, vE, psiN, wN, uN, vN)
    LInfErr[i,0] = psiInf
    LInfErr[i,1] = wInf
    LInfErr[i,2] = uInf
    LInfErr[i,3] = vInf
    #
    # Close the data file.
    fid.close()

# Calculate the error orders.
L2Order = np.zeros((nRest,4), dtype=np.float64)
LInfOrder = np.zeros((nRest,4), dtype=np.float64)
for i in range(1,nRest):
    dtDiff = np.log(dt[i-1] - dt[i])
    #
    # L2 errors.
    psiL2Diff = np.log(L2Err[i-1,0] - L2Err[i,0])
    wL2Diff = np.log(L2Err[i-1,1] - L2Err[i,1])
    uL2Diff = np.log(L2Err[i-1,2] - L2Err[i,2])
    vL2Diff = np.log(L2Err[i-1,3] - L2Err[i,3])
    L2Order[i,0] = psiL2Diff/dtDiff
    L2Order[i,1] = wL2Diff/dtDiff
    L2Order[i,2] = uL2Diff/dtDiff
    L2Order[i,3] = vL2Diff/dtDiff
    #
    # LInf errors
    psiLInfDiff = np.log(LInfErr[i-1,0] - LInfErr[i,0])
    wLInfDiff = np.log(LInfErr[i-1,1] - LInfErr[i,1])
    uLInfDiff = np.log(LInfErr[i-1,2] - LInfErr[i,2])
    vLInfDiff = np.log(LInfErr[i-1,3] - LInfErr[i,3])
    LInfOrder[i,0] = psiLInfDiff/dtDiff
    LInfOrder[i,1] = wLInfDiff/dtDiff
    LInfOrder[i,2] = uLInfDiff/dtDiff
    LInfOrder[i,3] = vLInfDiff/dtDiff

# Write out the error analysis to file.
hdr = '    dt     ' + \
      '  psiL2Err.  ' + \
      '  psiL2Ord.  ' + \
      '   wL2Err.   ' + \
      '   wL2Ord.   ' + \
      '   uL2Err.   ' + \
      '   uL2Ord.   ' + \
      '   vL2Err.   ' + \
      '   vL2Ord.   ' + \
      '  psiLIErr.  ' + \
      '  psiLIOrd.  ' + \
      '   wLIErr.   ' + \
      '   wLIOrd.   ' + \
      '   uLIErr.   ' + \
      '   uLIOrd.   ' + \
      '   vLIErr.   ' + \
      '   vLIOrd.   '
np.savetxt('Error.dat', np.column_stack((dt, \
                                         L2Err[:,0], L2Order[:,0], \
                                         L2Err[:,1], L2Order[:,1], \
                                         L2Err[:,2], L2Order[:,2], \
                                         L2Err[:,3], L2Order[:,3], \
                                         LInfErr[:,0], LInfOrder[:,0], \
                                         LInfErr[:,1], LInfOrder[:,1], \
                                         LInfErr[:,2], LInfOrder[:,2], \
                                         LInfErr[:,3], LInfOrder[:,3])), \
           header=hdr, fmt='%11.6e')
