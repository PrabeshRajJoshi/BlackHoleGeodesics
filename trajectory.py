# -*- coding: utf-8 -*-
"""
Created on Fri Apr 03 21:35:43 2015

@author: prabesh-sv
"""

import numpy as np
from numpy import sin,cos,pi,abs,sqrt,arcsin,isfinite
from math import isnan

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
#import matplotlib.cm as cmx

import matplotlib
matplotlib.rcParams.update({'axes.labelsize': 18})

## ========================================================================

"""
RK4             : carry out the Runga Kutta 4 method of integration for solution
scatter_plot    : plot a 3D scatter diagram showing the complete particle trajectory.
single_plot     : plot and save 3D diagrams, each with position of particle at specified instance
check_solutions : plot the coordinates solutions (t,r,theta,phi) as functions of proper time \tau
z-rho_projection: plot 2D projection in plane[z=r*cos(theta) vs \rho = r*sin(theta)], and z=0 plane
incoming_particle: find solution for particle going towards the BH (if 0, outgoing particle solution is found)
phase_diagram   : plot phase space diagram and effective potentials for 'r' and '\theta' coordinates.
"""
criteria = { 'RK4'              :1,
             'scatter_plot'     :1,
             'single_plot'      :0,
             'check_solutions'  :1,
             'z-rho_projection' :1,
             'phase_diagram'    :1}
## ========================================================================

## ======================== FUNCTION DEFINITIONS ==========================

## ------------------------------------------------------------------------
def sphr2cart(sphr):
    """
    Function converts spherical coordinates sphr= [r, theta, phi] and returns cartesian coordinates [x, y, z].
     r = radial coordinate.
     theta = co-latitude.
     phi = azimuthal angle.
    """
    r = np.array(sphr[0])
    th = np.array(sphr[1])
    ph = np.array(sphr[2])

    x=r*sin(th)*cos(ph)
    y=r*sin(th)*sin(ph)
    z=r*cos(th)

    return [x,y,z]
## ------------------------------------------------------------------------

## functions in the first order ODE for each coordinate
## ------------------------------------------------------------------------
def ftime(rf):
    """
    Function requires global variables- m: black hole mass paramter, q : black hole charge parameter,
    E : Total energy of particle
    """
    Sig = 1 - 2*m/rf + q**2/rf**2
    return E/Sig
## ------------------------------------------------------------------------
def fphi(rf,thf):
    """
    Function requires global variables- b: cosmic string deficit angle paramter,
    Lz : z component of angular momentum of particle
    """
    denom = rf**2 * b* (sin(thf))**2
    return Lz/denom
## ------------------------------------------------------------------------
def fr1(rf):
    """
    Function for increasing r, requires global variables- m: black hole mass paramter, q : black hole charge parameter,
    E : Total energy of particle, Q : Carter constant ,
    ep : constant in lagrangian, and hamilton's principal function
    """
    Sig = 1 - 2*m/rf + q**2/rf**2
    R2 = Sig*(ep-Q/rf**2) + E**2
    return sqrt(abs(R2))
## ------------------------------------------------------------------------
def fr2(rf):
    """
    Function for decreasing r, requires global variables- m: black hole mass paramter, q : black hole charge parameter,
    E : Total energy of particle, Q : Carter constant ,
    ep : constant in lagrangian, and hamilton's principal function
    """
    Sig = 1 - 2*m/rf + q**2/rf**2
    R2 = Sig*(ep-Q/rf**2) + E**2
    return -sqrt(abs(R2))
## ------------------------------------------------------------------------
def ftheta1(rf,thf):
    """
    Function for increasing theta, requires global variables- b: cosmic string deficit angle paramter,
    Lz : z component of angular momentum of particle, Q : Carter Constant.
    Function increases value of theta
    """
    Theta2 = Q-(Lz/sin(thf))**2
    return sqrt(abs(Theta2))/rf**2
## ------------------------------------------------------------------------
def ftheta2(rf,thf):
    """
    Function for decreasing theta, requires global variables- b: cosmic string deficit angle paramter,
    Lz : z component of angular momentum of particle, Q : Carter Constant.
    Function decreases value of theta
    """
    Theta2 = Q-(Lz/sin(thf))**2
    return -sqrt(abs(Theta2))/rf**2
## ========================================================================

## ========================================================================
## -----------------------Constants for particle motion -------------------
L = 7.      ## Total angular momentum of particle
Q = L*L     ## Carter constant Q=L^2
m = 1.      ## mass parameter of BH
q = 0.6     ## charge parameter of BH
b = 1.0     ## cosmic string parameter
E = 1.5     ## Total Energy of particle
Lz = 2.     ## Angular momentum of particle along z-axis 
ep = -1.    ## constant of Lagrangian
## ------------------------------------------------------------------------

## ------------------------constants for simulation -----------------------
[ct0,r0,th0,ph0] = [0., 10.0, pi/3 , pi/2]
[pt_min,pt_max] = [0,20]    ## range of proper time
dpt = 0.001                 ## steps of proper time
## ------------------------------------------------------------------------

## ------------------------event horizon and constraint--------------------
ehr = round(m+sqrt(m**2-q**2),4) ## event horizon radius
cyhr = round(m-sqrt(m**2-q**2),4)
print ("event horizon, ehr=%1.4f \ncauchy horizon, cyhr=%1.4f")%(ehr,cyhr)

k_inv = Lz/L ## for constraint on \theta from d\theta / d\phi calculation
print "Constraint for theta:  |sin(theta)|>=", k_inv

## ------------------------------------------------------------------------

## -----background shperical surface at event horizon / cauchy horizon ----
n = 100
u = np.linspace(0,2*pi,n)
v = np.linspace(0,pi,n)

## event horizon 1 (absolute horizon)
ehx = ehr * np.outer(cos(u), sin(v))
ehy = ehr * np.outer(sin(u), sin(v))
ehz = ehr * np.outer(np.ones(n), cos(v))

## event horizon 2 (Cauchy horizon)
ehx2 = cyhr * np.outer(cos(u), sin(v))
ehy2 = cyhr * np.outer(sin(u), sin(v))
ehz2 = cyhr * np.outer(np.ones(n), cos(v))
## ------------------------------------------------------------------------
## ========================================================================

## -------------------------- The phase diagrams --------------------------
if criteria['phase_diagram'] ==1:
    r_cp = np.linspace(.1,5.,10000)
    r_cp = np.round(r_cp,decimals=5)

    ## Array to store zeros of polynomial R = (dr/dtau)^2
    R_zeros = []

    ## for Reissner Nordstrom
    Sig_cp = 1 - 2*m/r_cp + q**2/r_cp**2
    R_cp = Sig_cp*(ep-Q/r_cp**2) + E**2
    sqrt_R = sqrt(R_cp)

    ## check for intervals in r when polynomial R remains positive
    for i in range(len(r_cp)-1):
        if isfinite(sqrt_R[i]) and isnan(sqrt_R[i+1]):
            print "Polynomial R, +ve to -ve at r = ",r_cp[i]
            R_zeros.append(r_cp[i])
        elif isnan(sqrt_R[i]) and isfinite(sqrt_R[i+1]):
            print "Polynomial R, -ve to +ve at r = ",r_cp[i]
            R_zeros.append(r_cp[i])

    ## To indicate the zeros of R =(dr/dtau)^2 in the phase diagram
    R_zeros = [round(val,4) for val in R_zeros]
    zeros_text = r"Zeros of $R, r_0 =$" +str(R_zeros).strip('[]')

    ## for comparison with polynomial of Schwarzschild BH
    Sig_cp_Sch = 1 - 2*m/r_cp
    R_cp_Sch = Sig_cp_Sch*(ep-Q/r_cp**2) + E**2

    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)
    ax.plot(r_cp,R_cp,'g',linewidth=3,label=r'$R$ for RN')
#    ax.plot(r_cp,sqrt_R,'g--',linewidth=3,label=r'$\sqrt{R}$ for RN')
#    ax.plot(r_cp,-sqrt_R,'r--',linewidth=3,label=r'-$\sqrt{R}$ for RN')
    ax.plot(r_cp,R_cp_Sch,'b',linewidth=3,label=r'$R$ for Schwarzschild')
#    ax.plot(r_cp,np.zeros(len(r_cp)))
    ax.set_xlabel(r'radial coordinate:$r$', color='blue')
    ax.set_ylabel(r'$R(r)$', color='blue')
    ax.legend()
    ax.grid()

    plt.suptitle("Parameters: $Q="+str(Q)+", m="+str(m)+", q="+str(q)+r", \beta="+str(b)+", E="+str(E)+", L_z="+str(Lz)+r", \epsilon="+str(ep)+"$",fontsize=18)
    plt.title(zeros_text)
    plt.show()
    plt.close()

## ========================================================================

## ------------------------ The effective potentials ----------------------
    Veff_r = 0.5*Sig_cp*(Q/r_cp**2 - ep)
    E_tot = 0.5* E**2 * np.ones(len(r_cp))

    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(121)
    ax.plot(r_cp,Veff_r,'b',linewidth=3,label=r'$V_{eff}$')
    ax.plot(r_cp,E_tot,'r',linewidth=3,label=r'$E_{tot}$')
    ax.plot(r_cp,R_cp,'g--',linewidth=3,label=r'$R$ for RN')
    ax.plot(r_cp,np.zeros(len(r_cp)))
    ax.set_xlabel(r'radial coordinate:$r$', color='blue')
    ax.set_ylabel(r'$V_{eff}(r)$', color='blue')
    ax.set_title(zeros_text)
    plt.legend()
    ax.grid()

    theta_cp = np.linspace(0.1,pi-0.1,1000)
    Veff_theta = (Lz/(sin(theta_cp)))**2
    Q_array = Q * np.ones(len(theta_cp))

    ax = fig.add_subplot(122)
    ax.plot(theta_cp,Veff_theta,'g',linewidth=3,label=r'$V_{eff}$')
    ax.plot(theta_cp,Q_array,'r',linewidth=3,label=r'$Q$')
#    ax.plot(theta_cp,np.zeros(len(theta_cp)))
    ax.set_xlabel(r'$\theta$', color='blue')
    ax.set_ylabel(r'$V_{eff}(\theta)$', color='blue')
    plt.legend()
    ax.grid()

    plt.suptitle("Parameters: $Q="+str(Q)+", m="+str(m)+", q="+str(q)+r", \beta="+str(b)+", E="+str(E)+", L_z="+str(Lz)+r", \epsilon="+str(ep)+"$",fontsize=18)
    plt.show()
    plt.close()

## ========================================================================

## --------------------------- The RK4 routine ----------------------------
if criteria['RK4']==1:
    ## solution array initialisation
    pt = [pt_min]
    ct = [ct0]
    r=[r0]
    th = [th0]
    ph=[ph0]
    nmax=int( (pt_max-pt_min)/dpt +1)
    print "Solution length should be :",nmax

    ## constraint on theta found from relation btw theta and phi
    inv_sin = arcsin(k_inv)
    ## constraint on r for +ve values of polynomial R (found from phase diagram)
    rmin = input("Enter minimum bound for r with +ve R: ")
    rmax = input("Enter maximum bound for r with +ve R: ")

    ## function to start evolution of theta with (needed due to constraint)
    ftheta = ftheta2
    ## function to start evolution of r with (needed due to constraint)
    fr = fr2

    ## RK4 algorithm
    for n in range(nmax):
        ptn = pt[n]
        ctn= ct[n]
        rn = r[n]
        thn = th[n]
        phn = ph[n]

        ct1 = dpt*ftime(rn)
        r1 = dpt*fr(rn)
        th1 = dpt*ftheta(rn,thn)
        ph1 = dpt*fphi(rn,thn)

        ct2 = dpt*ftime(rn+r1/2.0)
        r2 = dpt*fr(rn+r1/2.0)
        th2 = dpt*ftheta(rn+r1/2.0,thn+th1/2.0)
        ph2 = dpt*fphi(rn+r1/2.0, thn+th1/2.0)

        ct3 = dpt*ftime(rn+r2/2.0)
        r3 = dpt*fr(rn+r2/2.0)
        th3 = dpt*ftheta(rn+r2/2.0,thn+th2/2.0)
        ph3 = dpt*fphi(rn+r2/2.0, thn+th2/2.0)

        ct4 = dpt*ftime(rn+r3)
        r4 = dpt*fr(rn+r3)
        th4 = dpt*ftheta(rn+r3,thn+th3)
        ph4 = dpt*fphi(rn+r3, thn+th3)

        ct.append( ctn+(1.0/6.0)*(ct1+2.0*ct2+2.0*ct3+ct4) )
        r.append( rn+(1.0/6.0)*(r1+2.0*r2+2.0*r3+r4) )
        th.append( (thn+(1.0/6.0)*(th1+2.0*th2+2.0*th3+th4)) )
        ph.append( (phn+(1.0/6.0)*(ph1+2.0*ph2+2.0*ph3+ph4))%(2*pi) )
        pt.append(ptn+dpt)

        ## change theta function if solution goes out of constraint
        if not inv_sin < thn < pi-inv_sin:
            if thn < inv_sin:
                ftheta = ftheta1
            elif thn > pi-inv_sin:
                ftheta = ftheta2

        ## change r function if solution goes out of constraint
        if not rmin < rn < rmax:
            if rn < rmin:
                fr = fr1
            elif rn > rmax:
                fr = fr2

    ## change to cartesian coordinates
    [x,y,z] = sphr2cart([r,th,ph])
    
#    fig = plt.figure(figsize=(20,10))
#    ax = fig.add_subplot(111)
#    ax.plot(pt,x,'b.')
#    ax.plot(pt,y,'g--.')
#    ax.plot(pt,z,'r')
#    plt.show()
#    plt.close()
## ========================================================================

## --plot each coordinate solutions vs proper time affine parameter \tau --
if criteria['check_solutions']==1:    
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(221)
    ax.plot(pt,ct)
    ax.set_xlabel(r'affine parameter, $\tau$', color='blue')
    ax.set_ylabel(r'$t$', color='blue')

    ax = fig.add_subplot(222)
    ax.plot(pt,r)
    ax.plot(pt,np.ones(len(pt))*ehr,'b--',label='event horizon')
    ax.plot(pt,np.ones(len(pt))*cyhr,'g-.',label='cauchy horizon')
    ax.set_xlabel(r'$\tau$', color='blue')
    ax.set_ylabel(r'$r$', color='blue')
    ax.legend()

    ax = fig.add_subplot(223)
    ax.plot(pt,sin(th))
    ax.set_xlabel(r'$\tau$', color='blue')
    ax.set_ylabel(r'$\sin \theta$', color='blue')

    ax = fig.add_subplot(224)
    ax.plot(pt,sin(ph))
    ax.set_xlabel(r'$\tau$', color='blue')
    ax.set_ylabel(r'$\sin \phi$', color='blue')
    plt.show()
    plt.close()
## ========================================================================

## plot 2D projection in z=r*cos\theta and r*sin\theta axes, and z=0 plane 
if criteria['z-rho_projection']==1:
    fig = plt.figure(figsize=(20,10))

    ## projection in z=r*cos\theta and r*sin\theta axes
    ax = fig.add_subplot(121)
    x_cord = r*sin(th)
    ax.plot(x_cord,z)
    ax.set_xlabel(r'$\rho=r\sin(\theta)$', color='blue')
    ax.set_ylabel(r'$z=r\cos(\theta)$', color='blue')
    begin_x = x_cord[0]
    begin_z = z[0]
    ## mark initial position with green circle
    ax.plot([begin_x],[begin_z],'go')
    
    ## event horizon and cauchy horizon indicators
    th_val = np.linspace(0,pi,100)
    ehr_z = ehr*cos(th_val)
    ehr_x_axis = ehr*sin(th_val)
    cyhr_z = cyhr*cos(th_val)
    cyhr_x_axis = cyhr*sin(th_val)
    plt.fill(ehr_x_axis,ehr_z,'b',alpha=0.2)
    plt.fill(cyhr_x_axis,cyhr_z,'r',alpha=0.2)

    ## projection in x,y place
    ax = fig.add_subplot(122)
    ax.plot(x,y)
    ## mark initial position with green circle
    ax.plot([x[0]],[y[0]],'go')

    ## event horizon and cauchy horizon indicators
    phi_val = np.linspace(0,2*pi,100)
    ehr_x = ehr*cos(phi_val)
    ehr_y = ehr*sin(phi_val)
    cyhr_x = cyhr*cos(phi_val)
    cyhr_y = cyhr*sin(phi_val)
    plt.fill(ehr_x,ehr_y,'b',alpha=0.2)
    plt.fill(cyhr_x,cyhr_y,'r',alpha=0.2)

    ax.set_xlabel(r'$x$', color='blue')
    ax.set_ylabel(r'$y$', color='blue')

    plt.suptitle("Parameters: $Q="+str(Q)+", m="+str(m)+", q="+str(q)+r", \beta="+str(b)+", E="+str(E)+", L_z="+str(Lz)+r", \epsilon="+str(ep)+r", \tau_{final}="+str(pt[-1])+"$")
    plt.show()
    plt.close()
## ========================================================================

## ----------- 3D scatter plot of complete trajectory of particle ---------
if criteria['scatter_plot']==1:
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111, projection='3d')
    ## plot background shperical surface at event horizon
    ax.plot_surface(ehx, ehy, ehz,  rstride=4, cstride=4, linewidth=0,color='b', alpha=0.3)
    ax.plot_surface(ehx2, ehy2, ehz2,  rstride=4, cstride=4, linewidth=0,color='r', alpha=0.5)

    ## parameters to view the 3D plot
    #ax.view_init(elev=0., azim=0)

    ##scale the colormap with the values in time array
    ax.scatter(x, y,z, s=10,c=pt,cmap='jet',linewidth=0)
#    title_info = "Parameters: $Q="+str(Q)+", m="+str(m)+", q="+str(q)+", \beta="+str(b)+", E="+str(E)+", Lz="+str(Lz)+", \epsilon="+str(ep)+"$"
    ax.set_xlabel('x', color='blue')
    ax.set_ylabel('y', color='blue')
    ax.set_zlabel('z', color='blue')

    plt.title("Parameters: $Q="+str(Q)+", m="+str(m)+", q="+str(q)+r", \beta="+str(b)+", E="+str(E)+", L_z="+str(Lz)+r", \epsilon="+str(ep)+r", \tau_{final}="+str(pt[-1])+"$")
    plt.show()
    plt.close()
## ========================================================================

## ------------------ for individual trajectory instances -----------------
if criteria['single_plot']==1:
    #jet = cm = plt.get_cmap('jet')
    #cNorm  = colors.Normalize(vmin=ct[0], vmax=ct[-1])
    for i in range(len(x)):
        nfig = i%100
        if nfig==0:
            ##ax.plot([x[i]],[y[i]],[z[i]], markersize=80,c=ti[i], cmap='jet',norm=cNorm)
            #plt.show()
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(ehx, ehy, ehz,  rstride=4, cstride=4, linewidth=0,color='b', alpha=0.3)
            ax.plot([x[i]],[y[i]],[z[i]], markersize=10, markerfacecolor='g',marker='o', markeredgecolor='g')

            ## parameters to view the 3D plot
            ax.view_init(elev=30, azim=-60)
            # Set limits to the axes
            limit = 10
            ax.set_xlim(-limit, limit)
            ax.set_ylim(-limit, limit)
            ax.set_zlim(-limit, limit)

            ## Write the axis labels
            ax.set_xlabel('x', color='blue')
            ax.set_ylabel('y', color='blue')
            ax.set_zlabel('z', color='blue')
            plt.show()
            plt.savefig('plots/plot'+str(i/100)+'.png')
            plt.close()
## ========================================================================
## JOB DONE!!
print "Fear the Darkness!"