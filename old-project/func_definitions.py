import os
import matplotlib.pyplot as plt
from subprocess import call

## find the position of maximum and minimum potential in array V_eff##
def find_max_min(array):
    max_pot = max(array)
    min_pot = min(array)
    max_id = [i for i,val in enumerate(array) if abs(max_pot-val)<1e-6]
    min_id = [i for i,val in enumerate(array) if abs(min_pot-val)<1e-6]

    if (len(max_id) and len(min_id))==1:
        idx = max_id[0]
        r_vmax = r_array[idx]
        idx = min_id[0]
        r_vmin = r_array[idx]
        print "max. potential position: ",r_vmax
        print "min. potential position: ",r_vmin

## ---------------------------------------------------------------------- ##

## find the position of zero in the array effective potential ##
def find_zero(array):
    zero_id = [i for i,val in enumerate(array) if val <1e-18]
    r_zeros = [r_array[val] for val in zero_id]
    print "zero effective potential at position: ", r_zeros
## ---------------------------------------------------------------------- ##

## ---------------------Effective Potential Calculator ------------------ ##
def eff_pot_calculator(r_arr, Energy, Ang_Momentum, MH):

    ## some combination of constants used in calculations ##
    aE_sq = (a*Energy)**2
    aEL_sq = (a*Energy - Ang_Momentum)**2
    L_sq = Ang_Momentum**2
    ##----------------------------------------------------##

    for i in range(int(length)):
        r = r_arr[i]
        M = MH *(r**3/(r**3 + g**3))
        # M = MH  ## for Kerr-solution
        delta = r**2 - 2*M*r + a**2

        ef_pot = (1./r**2) * (L_sq - delta*epsilon - aE_sq - (2*M*aEL_sq)/r)
        # velocity[i] = Energy**2 - ef_pot
        eff_pot[i] = ef_pot
    return eff_pot
## ---------------------------------------------------------------------- ##

## function for plotting multiple dependent values
def make_plot1(x_arr,y1_arr, y2_arr, x_label,y_label, plot_title, count, group_plot_name):
    fig_name = str(count)+"_"+group_plot_name+".png"

    plt.figure(figsize=(12,7))

    plt.plot(r_array, y1_arr,'b', label='$E^2$')
    plt.plot(r_array, y2_arr,'g', label='$V_{eff}$')

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(plot_title)
    plt.legend(loc='best')
    plt.axis([0,max(r_array),-5,5])	## limits of the axes
    plt.savefig(fig_name)
    plt.close()

## ---------------------------------------------------------------------- ##

## function for normal x-y plotting
def make_plot2(x_arr,y_arr, x_label,y_label, plot_title, count, group_plot_name):
    fig_name = str(count)+"_"+group_plot_name+".png"

    plt.plot(x_arr, y_arr, 'g')
    plt.plot(x_arr,y_arr2,'g')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(plot_title)
    plt.axis([0,10,-10,10]) ## limits of the axes
    plt.grid()
    plt.savefig(fig_name)
    plt.clf()
## ---------------------------------------------------------------------- ##

def make_movie(group_plot_name):
    fig_name_end = "_"+group_plot_name+".png"
    movie_name = group_plot_name+".mp4"

    ## make video from images (to install ffmpeg: sudo apt-get install ffmpeg)
    print"Making movie(frame rate=5).\n"
    command_make_video = "ffmpeg -qscale 1 -r 5 -b 9600 -i %d"+fig_name_end+" "+movie_name
    call(command_make_video,shell=True)

    print "Removing images...\n"
    command = "rm -rf *"+fig_name_end
    call (command,shell=True)

def E_L_solution(r_arr):
    '''
    solution for energy and angular momentum are obtained from the following set of equations

    R:= r{dot}^2 = E^2 + (1/r^2) [delta*epsilon + (a*E)^2 -L^2 + (2M/r) (aE-L)^2] 

    Now, R = 0 and dR/dr = 0 give following two set of equations

    L^2 = delta*epsilon + E^2 (a^2+(3/5)r^2)   and    (5M/r) (aE-L)^2 - (E*r)^2 = 0

    First a quadratic equation of the form AX^2 + BX + C = 0 with X = E^2 
    Then, E=sqrt(X) is obtained and finally L can be obtained from E and r
    '''
    ## initialise the arrays
    E_arr = []
    L_arr = []

    for r in r_arr:
        ## calculate changing M and delta
        M = M_H *(r**3/(r**3 + g**3))
        # M = inf_m  ## for Kerr-solution
        delta = r**2 - 2*M*r + a**2

        ## intermediate term(s) used in the solution of E
        de = delta*epsilon
        C1 = 5.*M*de
        A1 = 10.*(a**2)*M + 3.*(r**2)*M - r**3

        ## Coefficient of X^2 = E^4
        A = A1**2 - 100.*(M**2)*(a**4) - 60.*(M**2)*(a**2)*(r**2)
        ## Coefficient of X = E^2
        B = 2.*A1*C1 - 100.*(M*a)**2*de
        ## Constant in quadratic equation
        C = C1**2
        ## discriminant for the quadratic equation for X with X=E^2
        discri = B**2 - 4*A*C

        if discri >= 0:
            ## calculate square root of discriminant
            C2 = discri**(0.5)
            ## calculate the roots
            X1 = (-B + C2)/(2.*A)
            X2 = (-B - C2)/(2.*A)
        else:
            print "Discriminant (B^2 - 4AC) is negative! Program exit from function :E_L_solution \n"
            break
            # sys.exit("Discriminant (B^2 - 4AC) is negative! Program exit from function :EL_solution")

        ## find the value of energy
        if X1>=0:
            E_val = X1**(0.5)
        elif X2>=0:
            E_val = X2**(0.5)
        else:
            print "Real value of energy is not found. Program exit from function :E_L_solution \n"
            break

        ## check the condition for angular momentum to be positive
        check_L = (E_val*a)**2 + (3./5.)*(E_val*r)**2
        if check_L > de:
            ## value of L can be calculated from value of E and r
            L_val = (de + check_L)**(0.5)

        else:
            print "Angular momentum is bound to be imaginary! Program exit from function :E_L_solution \n"
            break
            # sys.exit("Angular momentum is bound to be negative! Program exit from function :EL_solution")

        ## check for any NaN values before adding it to the energy and momentum array
        if math.isnan(E_val) == False:
            E_arr.append(E_val)
            L_arr.append(L_val)

    return [E_arr,L_arr]

## ---------------------------------------------------------------------- ##

##========================================================================##

