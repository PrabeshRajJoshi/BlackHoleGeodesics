# -*- coding: utf-8 -*-
"""
Prabesh Joshi, 06/2014
"""

import pickle
from sympy import symbols,diff,latex,solve,div,sin,sqrt,simplify,integrate,discriminant
from func_definitions import make_plot2, make_movie

## ================= SPECIFY COMPONENTS TO SOLVE FOR ================== ##
# INDICATE 0(do not calculate) OR 1(calculate) FOR THE FOLLOWING CATEGORIES:
# R_zeros : find zeros of the polynomial R(r)
# dR_zeros: find zeros of derivative dR/dr
# EL_sol :  find energy/angular momentum solutions
# triple_sol: find triple point solution for R = dR/dr = d^2R/dr^2 =0
# horizon:  calculate the horizons from the \Delta parameter
# EL_plot: execute the python program to plot the saved solutions
# rtp : find dr/dt and dr/d\phi solutions

solve_for = {'R_zeros'      :1,
             'dR_zeros'     :0,
             'EL_sol'       :0,
             'triple_sol'   :0,
             'horizon'      :0,
             'EL_plot'      :0,
             'rtp'          :0}
          
## ==================================================================== ##

## prefix for solution filenames ##
EL_sol_file_prefix = "EL_solutions"

## ==================================================================== ##


## Define the necessary symbols/variables in the function
## e: epsilon parameter, constant for lagrangian
## E,L: Energy and Angular Momentum along z-axis
## t,r,theta,phi : spherical coordinate variables t,r, theta and phi
## m : mass of particle
## q : charge of particle
## b : deficit angle parameter
E,L,t,r,theta,phi,e,m,q,b = symbols('E,L,t,r,theta,phi,e,m,q,b')

## additional dotted variables
## td = dt/ds, rd = dr/ds , thd = d theta/ds, phd = d phi/ds
td,rd,thd,phd = symbols('td,rd,thd,phd')

## Sigma parameter from metric terms g_tt and g_rr of RN BH
Sig = 1 - 2*m/r + q**2/r**2
#print latex(Sig)

## Lagrangian for RN BH
Lag = Sig*td**2 - rd**2/Sig - r**2*(thd**2 + b**2 * (sin(theta))**2 *phd**2)
#print latex(Lag)


### Hamilton Jacobi Approach
#Q,delta = symbols('Q,delta')
### Obtained Theta and R polynomials
#THETA = Q-(L/b*sin(theta))**2
#R = delta/Sig - Q/(Sig*r**2) - (E/Sig)**2
#print integrate(sqrt(THETA),theta)
##print latex(THETA)
##print latex(R)

#
### Polynomial R(r) for Regular BH
#R = E**2*r**2 + e*r**2 +2*M*(a*E-L)**2 * (r**2/(r**3+g**3))- 2*e*M*r**4/(r**3+g**3) + e*a**2 + E**2*a**2 - L**2
### transform R(r) to (g^3+r^3)R for easier differentiation.
#R = (r**3 + g**3)*(E**2*r**2 + e*r**2 + e*a**2 + E**2*a**2 - L**2) + 2*M*(a*E-L)**2 * r**2- 2*e*M*r**4
#
## Polynomial R(r) for Kerr BH: results in 3rd order polynomial
#R = E**2*r**2 + (r**2-2*M*r+a**2)*e + a**2*E**2 -L**2 + (2*M/r)*(a*E-L)**2
#
#print "\n                       The associated polynomial, R(r):\n", latex(R)
#
### find the first derivative dR/dr
#dR = diff(R,r)
#print "\n                       First derivative, dR:\n ", latex(dR)


## ======================= SOLVE FOR E & L ======================== ##
if solve_for['EL_sol'] == 1:
    ## solve R=0 and dR/dr =0 (circular orbits) for E and L
    EL_sol = solve([R,dR],[E,L])
    print "\n       Number of solution sets (E,L):\n ", len(EL_sol)

    ## save the solutions
    for i in range(len(EL_sol)):
        
        with open(EL_sol_file_prefix+str(i)+".dat","w") as ELfile:
            pickle.dump(EL_sol[i], ELfile)
        ELfile.close()
        
#        print "E"+str(i+1)+" = ", latex(EL_sol[i][0])
#        print "\n L"+str(i+1)+" = ", latex(EL_sol[i][1])
#        print "\n\n"
    
    #print "\n"
    #print EL_sol[0][0]+EL_sol[1][0]
    #print "\n\n"
    #print EL_sol[0][0]+EL_sol[2][0]
        
## ================================================================= ##

##==================== FIND ZEROS OF R AND dR ======================##
if solve_for['R_zeros'] == 1:
    R_zeros = solve(R,r)
    print "\n       zeros of polynomial R(r):\n", R_zeros
    
if solve_for['dR_zeros'] == 1:
    dR_zeros = solve(dR,r)     
    print "\n       zeros of polynomial dR(r):\n", dR_zeros
## ================================================================= ##

## ============== CALCULATE THE TRIPLE POINT FOR R(r) ============== ##

if solve_for['triple_sol'] == 1:
    ## find d^2R/dr^2
    d2R = diff(R,r,2)
    print "\n       Second derivative, d2R:\n", latex(d2R)

    ## solve d^2R/dr^2 = 0 for the triple point value, r=r(a,g)
    triple_sol = solve([R,dR,d2R],r)
    print "\n       Number of solutions: ", len(triple_sol)
    
    for i in range(len(triple_sol)):
        print "r"+str(i)+": ", latex(triple_sol[i])
        print "\n\n"

## ================================================================= ##

## ============== CALCULATE THE HORIZON FROM \Delta(r) ============== ##

if solve_for['horizon'] == 1:
    a=0.5
    g=0.2
    M=1
    Delta = (r**3+g**3)*(r**2 + a**2) - (2*M*r**4)
    horizon = solve(Delta,r)
    print "No. of horizon solutions:", len(horizon)
    for i in range(len(horizon)):
        print "r_"+str(i)+":  ",horizon[i]
    quit()

## ================================================================= ##


## ====================== MAKE E/L PLOTS =========================== ##
if solve_for['EL_plot'] == 1:
    
    ##----- variables to make simulation video for a changing variable -----##
    sim_var = "M_H"
    delta_sim_var = 0.5     ## step size to change the variable
    min_sim_var = 0.5   
    max_sim_var = 5.0       ## minimum and maximum value of changing variable

    execfile("plot_EL_solutions.py")
## ================================================================= ##

## ====================== FIND dr/dt AND dr/d\phi ========================== ##
if solve_for['rtp'] == 1:
    
    ## Use definition of dt/d\tau , d\phi/d\tau and dr/d\tau ##
    Delta = (r**3+g**3)*(r**2 + a**2) - (2*M*r**4)
    dtdtau = ((r**2 + a**2 + (2*a**2*M)/r)*E - (2*a*M/r)*L)/Delta
    dpdtau = ((1 - 2*M/r)*L + (2*a*M/r)*E)/Delta
    drdtau_sq = E**2 + ( Delta*e + a**2*E**2 - L**2 + (2*M/r)*(a*E-L)**2 )/r**2
    
#    print "\n           Delta \n",latex(dtdtau**2)
#    exit()
    
    ## Calculate dr/d\phi and dr/dt ##
    drdt_sq = drdtau_sq/dtdtau**2
    drdphi_sq = drdtau_sq/dpdtau**2
    
    print "\n           dr_dt \n",latex(drdt_sq)
    print "\n           dr_dphi \n",latex(drdphi_sq)
    
    execfile("plot_drdt_drdphi.py")
    
#    ## Calculate dr/d\phi and dr/dt in quotient and remainder format ##
#    drdt = div(drdtau,dtdtau)
#    drdphi = div(drdtau,dpdtau)
#    
#    print "\n           dr_dt \n",latex(drdt)
#    print "\n           dr_dphi \n",latex(drdphi)
    
## ================================================================= ##


print "Task completed!\n"