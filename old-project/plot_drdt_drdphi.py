# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 15:27:39 2014

@author: prabesh-sv
"""


from sympy import lambdify
import numpy as np
import pickle
from matplotlib import pyplot as plt

## ================== DEFINE PARAMETER VALUES ==================== ##
values = {M:1.0, a:1.0, g:1.0, e:-1.0, E:0.9 , L:0.9}
arr_len = 1000
r_min = 0.5
r_max = 10
r_arr = np.linspace(r_min,r_max,arr_len)
## ==================================================================== ##

## drdt_sq and drdphi_sq are defined in analytic_solutions.py ##

## enter values of the parameters ##
drdt2 = drdt_sq.subs(values)
drdphi2 = drdphi_sq.subs(values)

### check the arrays ##
#print "\n           dr_dt \n",drdt2
#print "\n           dr_dphi \n",drdphi2
#quit()

## Make drdt2 and drdphi2 as functions of r (so that numpy array r_arr can be used)##
drdt2_r = lambdify(r,drdt2,'numpy')
drdphi2_r = lambdify(r,drdphi2,'numpy')

## Find numerical values of two the two arrays ##
drdt2 = drdt2_r(r_arr)
drdphi2 = drdphi2_r(r_arr)

## Plot the solutions ##

plot_title = "$M_H="+str(values[M])+",\; a="+str(values[a])+",\; g="+str(values[g])+",\; \epsilon="+str(values[e])+",\; E="+str(values[E])+",\; L="+str(values[L])+"$"

plt.plot(r_arr,drdt2,'b-', label="$(dr/dt)^2$")
plt.plot(r_arr,drdphi2,'g-', label="$(dr/d\phi)^2$")
plt.xlabel("$r$")
plt.ylabel("$(dr/dt)^2$ or $(dr/d\phi)^2$")
plt.title(plot_title)
plt.grid()
plt.legend(loc="best")
plt.show()
plt.clf()

