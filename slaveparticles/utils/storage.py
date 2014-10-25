# -*- coding: utf-8 -*-
"""
@author: oscar
Created on Mon May 19 18:21:43 2014
netCDF storage
"""
import time
import numpy as np
from netCDF4 import Dataset
#rootgrp = Dataset('CF_3b5e.nc', 'w', format='NETCDF4')
#cf0 = rootgrp.createGroup('degenerate')

#rootgrp.close()
#rootgrp.description = '3 Band system, with one lifted. 5 electrons'
#rootgrp.history = 'Created ' + time.ctime(time.time())
#rootgrp.source = 'netCDF4 python module tutorial'
#rootgrp.close()

def setgroup(group):

    group.history = 'Created ' + time.ctime(time.time())
    group.createDimension('Uintera', None)
    group.createDimension("slaves", 6)
    group.createDimension('MFields', 2)

    uin = group.createVariable('Uintera', 'f8', ('Uintera',))
#    slave= group.createVariable('slaves','u1', ('slaves',))
    zet = group.createVariable('Quasiparticle', 'f8', ('Uintera', 'slaves',))
    lag = group.createVariable('LagMulti', 'f8', ('Uintera', 'slaves',))
    orb = group.createVariable('OrbitalEnergy', 'f8', ('Uintera', 'slaves',))
    nup = group.createVariable('Upper orb pop', 'f8', ('Uintera',))
    tcf = group.createVariable("Crystal Field", 'f8')
    mfl = group.createVariable('Mean Fields', 'f8', ('Uintera', 'MFields', 'slaves',))

    return uin, zet, lag, orb, nup, tcf, mfl

def storegroup(variables, u_span, zet, lam, mu, nup, cf, mean_f):

    suin, szet, slag, sorb, snup, stcf, smfl = variables

    suin[:] = u_span
    szet[:] = np.asarray(zet)
    slag[:] = np.asarray(lam)
    sorb[:] = np.asarray(mu)
    snup[:] = np.asarray(nup)
    stcf[:] = cf
    smfl[:] = mean_f

def readgroup(group):

    uin = group.variables['Uintera'][:]
    zet = group.variables['Quasiparticle'][:]
    lag = group.variables['LagMulti'][:]
    orb = group.variables['OrbitalEnergy'][:]
    nup = group.variables['Upper orb pop'][:]
    tcf = group.variables['Crystal Field'][:]
    mfl = group.variables['Mean Fields'][:]

    return uin, zet, lag, orb, nup, tcf, mfl