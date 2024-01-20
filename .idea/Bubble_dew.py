import math
from scipy.optimize import fsolve
import numpy as np
import random
def antoine_eqn(T, A, B, C):
    return math.e**(A - B / (T + C))
def temp_vapour_pressure (A,B,C,P):
    return (B/(A-np.log(P)))-C

def raoults_law(x, P_sat):
    return sum(x * P_sat)

def bubble_point_pressure(T, x, A, B, C):
    # Function to solve for bubble point pressure
    def func(P):
        P_sat = np.array([antoine_eqn(T, a, b, c) for a, b, c in zip(A, B, C)])
        return raoults_law(x, P_sat) - P

    # Initial guess for pressure
    P_guess = 100.0  # Initial guess in bar

    # Solve for bubble point pressure using fsolve
    P_bubble = fsolve(func, P_guess,xtol=1e-6)[0]
    P_sat = np.array([antoine_eqn(T, a, b, c) for a, b, c in zip(A, B, C)])
    yi= (x * P_sat)/P_bubble
    return print('The bubble pressure is :- ',P_bubble),print ('The Vapor phase composition is ',yi)

def dew_point_pressure(T, y, A, B, C):
    # Function to solve for dew point pressure
    def func(P):
        P_sat = np.array([antoine_eqn(T, a, b, c) for a, b, c in zip(A, B, C)]) #P_vap at a prescribed t
        return (np.sum(y/P_sat))**-1 - P # Function to find root of

    # Initial guess for pressure
    P_guess = 100.0  # Initial guess in bar

    # Solve for dew point pressure using fsolve
    P_dew = fsolve(func, P_guess,xtol=1e-10)[0]
    P_sat = np.array([antoine_eqn(T, a, b, c) for a, b, c in zip(A, B, C)])
    xi=(y*P_dew)/P_sat
    return print('The dew point pressure is :- ',P_dew),print ('The liquid phase composition is :- ',xi)


def bubble_point_temperature(P, x, A, B, C):
    # Function to solve for bubble point temperature
    def func(T):
        P_sat = np.array([antoine_eqn(T, a, b, c) for a, b, c in zip(A, B, C)])
        return sum(raoults_law(x, P_sat) - P)

    # Initial guess for temperature
    temperature_sat = temp_vapour_pressure(A,B,C,P)
    T_guess =np.sum(x*temperature_sat)   # Initial guess is based on the huristic that the bubble temperature of the mixture will be dominated by the majority liquid component

    # Solve for bubble point temperature using fsolve
    T_bubble = fsolve(func, T_guess,xtol=1e-6)[0] # solving the non-linear equation
    P_sat = np.array([antoine_eqn(T_bubble, a, b, c) for a, b, c in zip(A, B, C)])
    yi=(x*P_sat)/P
    return print('The Bubble point temperature is ',T_bubble),print('Vapour feed composition:- ',yi)

def dew_point_temperature(P, y, A, B, C):
    # Function to solve for dew point temperature
    def func(T):
        P_sat = np.array([antoine_eqn(T, a, b, c) for a, b, c in zip(A, B, C)])
        z=y/P_sat
        return (np.sum(z))**-1 - P

    # Initial guess for temperature
    temperature_sat = temp_vapour_pressure(A, B, C, P)
    T_guess = np.sum(y * temperature_sat)# initial guess

    # Solve for dew point temperature using fsolve
    T_dew = fsolve(func, T_guess,xtol=1e-6)[0]
    P_sat = np.array([antoine_eqn(T_dew, a, b, c) for a, b, c in zip(A, B, C)])
    xi=(y*P)/P_sat
    return print('The Dew point temperature is :- ',T_dew), print('\nThe liquid phase composition is :- ',xi)
# Example usage:
T_specified = 300.0  # Specified temperature in kelvin
x_mixture = np.array([0.4, 0.3,0.3])  # Mole fractions of components in liquid phase
y_mixture = np.array([0.4, 0.3,0.3])  # Mole fractions of components in vapor phase
A_constants = [9.2806,9.3935,9.5188]  # Antoine equation constants for component 1 and 2
B_constants = [2788.51,3096.52,3366.99]
C_constants = [-52.36,-53.67,-58.04]

P_bubble_point = bubble_point_pressure(T_specified, x_mixture, A_constants, B_constants, C_constants)
P_dew_point = dew_point_pressure(T_specified, y_mixture, A_constants, B_constants, C_constants)
dew_point_temperature(0.029032709542123065,y_mixture,A_constants,B_constants,C_constants)

