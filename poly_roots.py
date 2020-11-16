# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 19:45:13 2015

@author: prabesh-sv
"""


from sympy import symbols,diff,latex,solve,div,sin,sqrt,simplify,integrate,discriminant
import pickle

E,L,t,r,e,m,q,b,Q = symbols('E,L,t,r,e,m,q,b,Q')

## prefix for solution filenames ##
EL_sol_file_prefix = "EL_solutions"

solve_for = {'R_zeros'      :0,
             'dR_zeros'     :0,
             'EL_sol'       :1,
             'triple_sol'   :0,
             'horizon'      :0,
             'EL_plot'      :1,
             'rtp'          :0}
             
             
Sig = 1 - 2*m/r + q**2/r**2
R2 = Sig*(e-Q/r**2) + E**2  ## R^2 := (dr/d\tau)^2


if solve_for['R_zeros'] == 1:
    R_zeros = solve(R2,r)
    print("\nZeros of polynomial R(r):\n", R_zeros)

## ======================= SOLVE FOR E & L ======================== ##
if solve_for['EL_sol'] == 1:
    ## solve R=0 and dR/dr =0 (circular orbits) for E and L
    
    dR2 = diff(R2,r)
    EL_sol = solve([R2,dR2],[E,L])
    print( "\nNumber of solution sets (E,L): %i \n" %(len(EL_sol)) )

    ## save the solutions
    for i in range(len(EL_sol)):
        # open the writeable file in binary mode for pickle to work!
        with open(EL_sol_file_prefix+str(i)+".dat","wb") as ELfile:
            pickle.dump(EL_sol[i], ELfile)
        ELfile.close()
        
#        print "E"+str(i+1)+" = ", latex(EL_sol[i][0])
#        print "\n L"+str(i+1)+" = ", latex(EL_sol[i][1])
#        print "\n\n"
        print(EL_sol[i])
    
    #print "\n"
    #print EL_sol[0][0]+EL_sol[1][0]
    #print "\n\n"
    #print EL_sol[0][0]+EL_sol[2][0]
        
## ================================================================= ##
        
if solve_for['EL_plot']==1:
    
    from sympy import lambdify
    import numpy as np
    from matplotlib import pyplot as plt
    
    ## ================== DEFINE PARAMETER VALUES ==================== ##
    values = {m:1.0, q:1.0, Q:1.0, e:-1.0}
    arr_len = 1000
    r_min = 0.5
    r_max = 5
    r_arr = np.linspace(r_min,r_max,arr_len)
    ## ==================================================================== ##
    
    for i in range(len(EL_sol)):
        with open(EL_sol_file_prefix+str(i)+".dat","rb") as ELin:
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
    
    
    plot_title = "$m="+str(values[m])+",\; q="+str(values[q])+",\; Q="+str(values[Q])+",\; \epsilon="+str(values[e])+"$"
    
    plt.plot(r_arr,Er0,'b-', label="$E_0$")
    plt.plot(r_arr,Er1,'g-', label="$E_1$")
    plt.xlabel("$r$")
    plt.ylabel("$E(r)$")
    plt.title(plot_title)
    plt.grid()
    plt.legend()
    plt.savefig("plots/E0_E1_r.png")
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

    plt.plot(Er0,Lr0,'b-', label="$E_0/L_0$")
    plt.plot(Er1,Lr1,'g-', label="$E_1/L_1$")
    plt.xlabel("$E(r)$")
    plt.ylabel("$L(r)$")
    plt.title(plot_title)
    plt.grid()
    plt.legend()
    plt.savefig("plots/phase_EL1.png")
    plt.clf()
    
    '''
    plt.plot(r_arr,Er2,'b-', label="$E_2$")
    plt.plot(r_arr,Er3,'g-', label="$E_3$")
    plt.xlabel("$r$")
    plt.ylabel("$E(r)$")
    plt.title(plot_title)
    plt.grid()
    plt.legend()
    plt.savefig("plots/E2_E3_r.png")
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
    
    
    plt.plot(Er2,Lr2,'b-', label="$E_2/L_2$")
    plt.plot(Er3,Lr3,'g-', label="$E_3/L_3$")
    plt.xlabel("$E(r)$")
    plt.ylabel("$L(r)$")
    plt.title(plot_title)
    plt.grid()
    plt.legend()
    plt.savefig("plots/phase_EL2.png")
    plt.clf()
    '''
    
    ## END ##