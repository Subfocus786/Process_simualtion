import math
from scipy.optimize import fsolve
import numpy as np
import random
from Flash_code import flash_case3
def antoine_eqn(T, A, B, C):
    return 10**(A - B / (T + C))
def antonie_eq_e (T, A, B, C):
    return math.e**(A - B / (T + C))
def temp_vapour_pressure (A,B,C,P):
    return (B/(A-math.log(P)))-C

def raoults_law(x, P_sat):
    return sum(x * P_sat)

def bubble_point_pressure(T, x, A, B, C):
    # Function to solve for bubble point pressure
    def func_P(P):
        P_sat = antoine_eqn(T,A,B,C)

        return raoults_law(x, P_sat) - P

    # Initial guess for pressure
    P_guess = random.uniform(a=1.0,b=500.0)  # Initial guess in bar

    # Solve for bubble point pressure using fsolve
    P_bubble = fsolve(func_P, P_guess,xtol=1e-6)[0]
    P_sat = antoine_eqn(T,A,B,C)
    yi= (x * P_sat)/P_bubble
    return P_bubble #print ('The Vapor phase composition is ',yi)

def dew_point_pressure(T, y, A, B, C):
    # Function to solve for dew point pressure
    def func(P):
        P_sat = np.array([antoine_eqn(T, a, b, c) for a, b, c in zip(A, B, C)]) #P_vap at a prescribed t
        return (np.sum(y/P_sat))**-1 - P # Function to find root of

    # Initial guess for pressure
    P_guess = random.uniform(1.0,1000.0)  # Initial guess in bar

    # Solve for dew point pressure using fsolve
    P_dew = fsolve(func, P_guess,xtol=1e-10)[0]
    P_sat = np.array([antoine_eqn(T, a, b, c) for a, b, c in zip(A, B, C)])
    xi=(y*P_dew)/P_sat
    return P_dew #print ('The liquid phase composition is :- ',xi)


def bubble_point_temperature(P, x, A, B, C):
    # Function to solve for bubble point temperature
    def func(T):
        P_sat = antonie_eq_e(T,A,B,C)
        return raoults_law(x, P_sat) - P

    # Initial guess for temperature
    temperature_sat = np.array( [ temp_vapour_pressure(A=a,B=b,C=c,P=P) for a, b, c in zip(A, B, C)])
    #print ('this is the saturated temp array',temperature_sat)
    T_guess =np.sum(x*temperature_sat)   # Initial guess is based on the huristic that the bubble temperature of the mixture will be dominated by the majority liquid component
    np.array(T_guess)
    #print(T_guess)
    # Solve for bubble point temperature using fsolve
    T_bubble = fsolve(func,T_guess,xtol=1e-6)[0] # solving the non-linear equation
    P_sat = antoine_eqn(T_bubble,A,B,C)
    yi=(x*P_sat)/P
    return T_bubble #print('\nThe Bubble point temperature is ',T_bubble),print('\n Vapour feed composition:- ',yi)

def dew_point_temperature(P, y, A, B, C):
    # Function to solve for dew point temperature
    def func(T):
        P_sat = antoine_eqn(T,A,B,C)
        z=y/P_sat
        return (np.sum(z))**-1 - P

    # Initial guess for temperature
    temperature_sat = random.uniform(50.0,400.0)
    T_guess = np.sum(y * temperature_sat)# initial guess

    # Solve for dew point temperature using fsolve
    T_dew = fsolve(func, T_guess,xtol=1e-7)[0]
    P_sat = np.array([antoine_eqn(T_dew, a, b, c) for a, b, c in zip(A, B, C)])
    xi=(y*P)/P_sat
    return print('\nThe Dew point temperature is :- ',T_dew), print('\nThe liquid phase composition is :- ',xi)


'''Note :- always check antonie eq form'''
"""A=[15.9008,16.0137,16.1156]
Ai=np.array(A)
B =[2788.51,3096.52,3395.57]
Bi=np.array(B)
C =[-52.34,-53.67,-59.44]
Ci=np.array(C)
fi=np.array([30,50,40])
dew_point_temperature(750,fi,Ai,Bi,Ci)"""