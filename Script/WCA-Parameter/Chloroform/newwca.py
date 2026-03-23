import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
from scipy import optimize
from scipy.optimize import brentq
import sympy
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as uvs

import warnings

# Suppress specific UserWarning about numexpr and bottleneck versions
warnings.filterwarnings("ignore", message="Pandas requires version '2.8.4' or newer of 'numexpr'")
warnings.filterwarnings("ignore", message="Pandas requires version '1.3.6' or newer of 'bottleneck'")

# Suppress specific RuntimeWarning for division by zero
warnings.filterwarnings("ignore", message="divide by zero encountered in divide")


tempval = sys.argv[1]
temp = float(tempval)
ktemp = 0.5922/298.0*temp
beta = 1/ktemp
boxlength = float(sys.argv[2])
rho_HS = 1000/(boxlength**3)

xval=[]
yval=[]
#f =open("./fitting_MeOH_MeOH.table",'r')
fname = str(tempval)+".out"
f =open(fname,'r')
i=0
for line in f:
    if i > 0:
        line_e = line.split()
        xval.append(float(line_e[1]))
        yval.append(float(line_e[2]))
    i+=1
spl = uvs(xval,yval,k=5)

xval = np.array(xval)
yval = np.array(yval)
x_range = np.linspace(min(xval), max(xval), 1000)
y_spl = spl(x_range)
min_index = np.argmin(y_spl)
min_x = x_range[min_index]
min_y = y_spl[min_index]

# Shift the interaction so that the minimum becomes zero
yval_shifted = yval - min_y
spl_shifted = uvs(xval, yval_shifted, k=5)

def get_cavity_function(r, d, eta, verbose=False):
    # Return zeroth and first shell of cavity function
    from numpy import exp, cos, sin
    x = r/d

    # Zeroth shell
    lambda_1 = (1+2*eta)**2/(1-eta)**4
    lambda_2 = -(2+eta)**2/(4*(1-eta)**4)
    if verbose:
        print('    Zeroth shell parameters')
        print(f'lambda_1 = {lambda_1:0.4f}')
        print(f'lambda_2 = {lambda_2:0.4f}')
    x = r/d
    c = (-lambda_1-6*eta*lambda_2*x-0.5*eta*lambda_1*x**3)*np.heaviside(1-x, 0.5)
    y = -c

    # First shell
    f = 3 + 3*eta-eta**2
    arg = (2*eta**4/f**2+1)**0.5
    y_plus = (2*eta*f)**(1/3)*(arg+1)**(1/3)
    y_minus = (2*eta*f)**(1/3)*(arg-1)**(1/3)
    z_d = y_plus - y_minus
    z_s = y_plus + y_minus
    A = (-2*eta+z_d)/(1-eta)
    B = (-2*eta-0.5*z_d)/(1-eta)
    C = 3**(1/2)*z_s/(2*(1-eta))
    a_1 = -2*eta*(1-eta-3*eta**2) + (1-3*eta-4*eta**2)*z_d + (1+eta/2)*z_d**2
    a_1 = a_1/(3*(2*eta**2+z_d**2)*(1-eta)**2)
    a_2 = eta*(2+4*eta-3*eta**2)-(1-3*eta-4*eta**2)*z_d+2*(1+eta/2)*z_d**2
    a_2 = a_2/(3*(2*eta**2+z_d**2)*(1-eta)**2)
    a_3 = (1-3*eta-4*eta**2)*(4*eta**2+z_d**2)+eta*(2-5*eta**2)*z_d
    a_3 = a_3/(3**(1/2)*z_s*(2*eta**2+z_d**2)*(1-eta)**2)
    if verbose:
        print('    First shell parameters')
        print(f'A = {A:0.5f}')
        print(f'B = {B:0.5f}')
        print(f'C = {C:0.5f}')
        print(f'a_1 = {a_1:0.5f}')
        print(f'a_2 = {a_2:0.5f}')
        print(f'a_3 = {a_3:0.5f}')
    Chang_H1 = (a_1*exp(A*(x-1))+a_2*exp(B*(x-1))*cos(C*(x-1))+a_3*exp(B*(x-1))*sin(C*(x-1)))
    y += (Chang_H1/x)*(1-np.heaviside(1-x, 0.5))
    
    # Second shell
    f = 3 + 3 * eta - eta ** 2
    arg = (2 * eta ** 4 / f ** 2 + 1) ** 0.5
    y_plus = (2 * eta * f) ** (1 / 3) * (arg + 1) ** (1 / 3)
    y_minus = (2 * eta * f) ** (1 / 3) * (arg - 1) ** (1 / 3)
    z_d = y_plus - y_minus
    z_s = y_plus + y_minus
    A = (-2 * eta + z_d) / (1 - eta)
    B = (-2 * eta - 0.5 * z_d) / (1 - eta)
    C = 3 ** (1 / 2) * z_s / (2 * (1 - eta))
    a_1 = -2 * eta * (1 - eta - 3 * eta ** 2) + (1 - 3 * eta - 4 * eta ** 2) * z_d + (1 + eta / 2) * z_d ** 2
    a_1 /= (3 * (2 * eta ** 2 + z_d ** 2) * (1 - eta) ** 2)
    a_2 = eta * (2 + 4 * eta - 3 * eta ** 2) - (1 - 3 * eta - 4 * eta ** 2) * z_d + 2 * (1 + eta / 2) * z_d ** 2
    a_2 /= (3 * (2 * eta ** 2 + z_d ** 2) * (1 - eta) ** 2)
    a_3 = (1 - 3 * eta - 4 * eta ** 2) * (4 * eta ** 2 + z_d ** 2) + eta * (2 - 5 * eta ** 2) * z_d
    a_3 /= (3 ** (1 / 2) * z_s * (2 * eta ** 2 + z_d ** 2) * (1 - eta) ** 2)
    b_1 = -4/3*eta*((2*eta**2-z_d**2)*(1-6*eta-3*eta**2+20*eta**3+15*eta**4)+z_d*eta*(16+24*eta-21*eta**2-13*eta**3+21*eta**4))/((2*eta**2+z_d**2)**3*(1-eta)**2)
    b_2 = -1.0*b_1
    b_3 = 8*(3**(1/2))*eta**2*((4*eta**2+z_d**2)*(2-10*eta-24*eta**2+30*eta**3+79*eta**4+21*eta**5-17*eta**6)+2*z_d*eta*(16+40*eta-eta**2-50*eta**3+11*eta**4+52*eta**5+13*eta**6))/(z_s**3*(2*eta**2+z_d**2)**3*(1-eta)**2)
    b_4 = -2/3*eta*(2*eta*(10+28*eta+21*eta**2-13*eta**3-19*eta**4)+z_d*(2-12*eta-18*eta**2+28*eta**3+27*eta**4)+z_d**2*(4-6*eta-18*eta**2-7*eta**3))/((2*eta**2+z_d**2)**2*(1-eta)**3)
    b_5 = -4/3*eta**2*(4*(6-30*eta-82*eta**2+58*eta**3+222*eta**4+94*eta**5-25*eta**6)-z_d*(24-10*eta-164*eta**2-156*eta**3+22*eta**4+41*eta**5)-z_d**2*(10+32*eta+15*eta**2-31*eta**3-26*eta**4))/(z_s**2*(2*eta**2+z_d**2)**2*(1-eta)**3)
    b_6 = 4/(3)**(1/2)*eta**2*(24-10*eta-164*eta**2-156*eta**3+22*eta**4+41*eta**5-z_d*(10+32*eta+15*eta**2-31*eta**3-26*eta**4))/(z_s*(2*eta**2+z_d**2)**2*(1-eta)**3)
    Chang_01 =  b_1+b_2
    Chang_H2 = (b_1*exp(A*(x-2))+b_2*exp(B*(x-2))*cos(C*(x-2))+b_3*exp(B*(x-2))*sin(C*(x-2))+b_4*(x-2)*exp(A*(x-2))+b_5*(x-2)*exp(B*(x-2))*cos(C*(x-2))+b_6*(x-2)*exp(B*(x-2))*sin(C*(x-2)))
    #correction = (Chang_0-Chang_01)*0.5
    y += (Chang_H2/x)*(1 - np.heaviside(2-x, 0.5))
    return y

def awc_integral(d, beta, verbose=False):
    r = np.linspace(xval[0],20,1000)
    #if 16.0 > 3*d:
    #    r = np.linspace(xval[0],16.0,100) #2**18)
    #else:
    #    r = np.linspace(xval[0],3*d,1000)
            
    #r = x*2**(1/6)
    #d = f_d*2**(1/6)
    eta = np.pi/6*rho_HS*d**3
    print_st = "d:%e eta:%e \n" %(d, eta)
    #print(print_st)
    y_d = get_cavity_function(r, d, eta)

    v_d = 1/np.heaviside(r-d, 0.5)-1
    #v_d = np.zeros_like(r)
    # Use np.where to handle different cases
    #v_d = np.where(r > d, 1/np.heaviside(r-d, 0.5) - 1, float('inf'))
    # Optionally, handle the case when r == d separately if needed
    #v_d = np.where(r == d, 1, v_d)

    e_d = np.exp(-beta*v_d)
    v_0 = spl_shifted(r)*np.heaviside(min_x-r,0.5) # (4*r**-12-4*r**-6 + 1)*np.heaviside(2**(1/6)-r, 0.5)
    e_0 = np.exp(-beta*v_0)
    Delta_e = e_0-e_d
    if verbose:
        plt.figure()
        plt.title(f'beta = {beta}, T={1/beta}')
        plt.plot(r, y_d)
        plt.plot(r, Delta_e)
        plt.plot(r, y_d*Delta_e*r*r)
        #plt.plot([2**(1/6), 2**(1/6)], [0, 5])
        plt.xlabel('Pair distance, $r$')
        plt.legend(['y_d', 'Delta_e', 'y_d*Delta_e', r'$2^{1/6}\simeq1.12$'])
        plt.show()
    return np.trapz(y_d*Delta_e*r*r, r)

def safe_awc_integral(f,beta):
    result = awc_integral(f,beta)
    if np.isnan(result):
        return 1e10
    else:
        if result > 1e10:
            return 1e10
        else:
            return result

fs = np.linspace(4, 6.0, 100) # d space 
awc = [awc_integral(f, beta) for f in fs]

awc = [safe_awc_integral(f, beta) for f in fs]
root = brentq(lambda f: safe_awc_integral(f,beta), 4, 6.0)
print('d where integral vanish:', root)
