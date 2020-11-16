# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 16:07:53 2014

@author: prabesh-sv
"""    

from sympy import lambdify
import numpy as np
import pickle
from matplotlib import pyplot as plt

## ================== DEFINE PARAMETER VALUES ==================== ##
values = {M:1.0, a:1.0, g:1.0, e:-1.0}
arr_len = 1000
r_min = 0.5
r_max = 5
r_arr = np.linspace(r_min,r_max,arr_len)
## ==================================================================== ##

#for i in range(4):
#    with open(EL_sol_file_prefix+str(i)+".dat","r") as ELin:
#        EL_sol = pickle.load(ELin)
#    ELin.close()
#    
##    print "\n                          Energy "+str(i)+"\n",latex(EL_sol[0])
##    print "\n                          A. Momentum "+str(i)+"\n",latex(EL_sol[1])
##    continue
#    
#    Er = EL_sol[0].subs(values)
#    Lr = EL_sol[1].subs(values)
#
#    Er1 = lambdify(r,Er,'numpy')
#    Lr1 = lambdify(r,Lr,'numpy')
#
#    Er = Er1(r_arr)
#    Lr = Lr1(r_arr)
##    print "\n 		Energy solution array \n ", Er
##    print "\n 		Angular Momentum solution array \n ", Lr
#    
##    Er = Er.subs({r:0.7})
##    Lr = Lr.subs({r:0.7})
##    print Er
##    print Lr
##    print "\n"
##    continue
#
#    plot_title = "$M_H="+str(values[M])+",\; a="+str(values[a])+",\; g="+str(values[g])+",\; \epsilon="+str(values[e])+"$"
#    
#    plt.plot(r_arr,Er)
#    plt.xlabel("$r$")
#    plt.ylabel("$E_"+str(i)+"(r)$")
#    plt.title(plot_title)
#    plt.grid()
#    plt.savefig("plots/E"+str(i)+"_r.png")
#    plt.clf()
#    
#    plt.plot(r_arr,Lr)
#    plt.xlabel("$r$")
#    plt.ylabel("$L_"+str(i)+"(r)$")
#    plt.title(plot_title)
#    plt.grid()
#    plt.savefig("plots/L"+str(i)+"_r.png")
#    plt.clf()
#    
#    plt.plot(Er,Lr)
#    plt.xlabel("$E_"+str(i)+"$")
#    plt.ylabel("$L_"+str(i)+"$")
#    plt.title(plot_title)
#    plt.grid()
#    plt.savefig("plots/E_L_plot"+str(i)+".png")
#    plt.clf()

for i in range(4):
    with open(EL_sol_file_prefix+str(i)+".dat","r") as ELin:
        EL_sol = pickle.load(ELin)
    ELin.close()
    
#    print "\n                          Energy "+str(i)+"\n",latex(EL_sol[0])
#    print "\n                          A. Momentum "+str(i)+"\n",latex(EL_sol[1])
#    continue
    
    E_sol = EL_sol[0].subs(values)
    L_sol = EL_sol[1].subs(values)

    Er = lambdify(r,E_sol,'numpy')
    Lr = lambdify(r,L_sol,'numpy')
    
    if i==0:
        Er0 = Er(r_arr)
        Lr0 = Lr(r_arr)
    elif i==1:
        Er1 = Er(r_arr)
        Lr1 = Lr(r_arr)
    elif i==2:
        Er2 = Er(r_arr)
        Lr2 = Lr(r_arr) 
    elif i==3:
        Er3 = Er(r_arr)
        Lr3 = Lr(r_arr)


plot_title = "$M_H="+str(values[M])+",\; a="+str(values[a])+",\; g="+str(values[g])+",\; \epsilon="+str(values[e])+"$"

plt.plot(r_arr,Er0,'b-', label="$E_0$")
plt.plot(r_arr,Er1,'g-', label="$E_1$")
plt.xlabel("$r$")
plt.ylabel("$E(r)$")
plt.title(plot_title)
plt.grid()
plt.legend()
plt.savefig("plots/E0_E1_r.png")
plt.clf()

plt.plot(r_arr,Er2,'b-', label="$E_2$")
plt.plot(r_arr,Er3,'g-', label="$E_3$")
plt.xlabel("$r$")
plt.ylabel("$E(r)$")
plt.title(plot_title)
plt.grid()
plt.legend()
plt.savefig("plots/E2_E3_r.png")
plt.clf()

plt.plot(r_arr,Lr0,'b-', label="$L_0$")
plt.plot(r_arr,Lr1,'g-', label="$L_1$")
plt.xlabel("$r$")
plt.ylabel("$L(r)$")
plt.title(plot_title)
plt.grid()
plt.legend()
plt.savefig("plots/L0_L1_r.png")
plt.clf()

plt.plot(r_arr,Lr2,'b-', label="$L_2$")
plt.plot(r_arr,Lr3,'g-', label="$L_3$")
plt.xlabel("$r$")
plt.ylabel("$L(r)$")
plt.title(plot_title)
plt.grid()
plt.legend()
plt.savefig("plots/L2_L3_r.png")
plt.clf()

plt.plot(Er0,Lr0,'b-', label="$E_0/L_0$")
plt.plot(Er1,Lr1,'g-', label="$E_1/L_1$")
plt.xlabel("$E(r)$")
plt.ylabel("$L(r)$")
plt.title(plot_title)
plt.grid()
plt.legend()
plt.savefig("plots/phase_EL1.png")
plt.clf()

plt.plot(Er2,Lr2,'b-', label="$E_2/L_2$")
plt.plot(Er3,Lr3,'g-', label="$E_3/L_3$")
plt.xlabel("$E(r)$")
plt.ylabel("$L(r)$")
plt.title(plot_title)
plt.grid()
plt.legend()
plt.savefig("plots/phase_EL2.png")
plt.clf()

## END ##