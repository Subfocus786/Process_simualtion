import numpy as np
import math
from scipy.optimize import fsolve

def vapour_pressure_temp (A,B,C,T):
    return math.e**(A - (B / (T + C))) #Antonie eq in 10 form

def eq1_case_P_guess(P,temp,A,B,C,fk,key,eta_n): #Sub-funtion for Case when Temp is specified and P is to guessed
    P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=temp) #Vapour pressure array from Antonie Eq
    # partition coefficent
    Kn = P_vap / P
    # calculation of the relative volatilities
    alpha_k_n = Kn / Kn[key - 1]
    # overhead split fraction calculation for all components
    Eta_k = (alpha_k_n * eta_n) / (1 + ((alpha_k_n - 1) * eta_n))
    # Molar fractional and flow calculation
    vk = Eta_k * fk
    lk = (1 - Eta_k) * fk
    Vtotal = np.sum(vk)
    Ltotal = np.sum(lk)
    xk = lk / Ltotal
    #yk = vk / Vtotal
    #print(P_vap);
    temp_pcalc = (P_vap * xk)
    p_calc = np.sum(temp_pcalc)
    #print(P,p_calc,temp,P_vap,xk)
    return (p_calc / P) - 1 #Total pressure calculation using raoult's law should be same as guess pressure
def eq1_case_t_guess(temp,P,A,B,C,fk,key,eta_n): #Sub-funtion for Case when P is specified and Temp is to guessed
    P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=temp) #Vapour pressure array from Antonie Eq at a guess temp
    Kn = P_vap / P
    # calculation of the relative volatilities
    alpha_k_n = Kn / Kn[key - 1]
    # overhead split fraction calculation
    Eta_k = (alpha_k_n * eta_n) / (1 + ((alpha_k_n - 1) * eta_n))
    # Steam table calculation
    vk = Eta_k * fk
    lk = (1 - Eta_k) * fk
    Vtotal = np.sum(vk)
    Ltotal = np.sum(lk)
    xk = lk / Ltotal
    yk = vk / Vtotal
    p_calc = np.sum(P_vap * xk)
    #print(P, p_calc, temp, P_vap)
    return (p_calc / P) - 1 #Total pressure calculation using raoult's law should be same as guess pressure

def flash_case1(A,B,C,fk,x=False): # Flash when Eta and P/T is given
    eta_n = float(input('Enter the overhead split fraction of the key component:  '))
    key = int(input("which component is the key"))
    if x==True:
        Temp_spec= float(input("Enter the specified temperature in Kelvin: "))
        P_guess1 = float(input("Enter the guess Pressure in mmHg"))
        result= fsolve(eq1_case_P_guess,P_guess1,args=(Temp_spec,A,B,C,fk,key,eta_n),xtol=1e-7) # using fsolves in scipy
        P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=Temp_spec)
        Kn = P_vap / result
        # calculation of the relative volatilities
        alpha_k_n = Kn / Kn[key - 1]
        # overhead split fraction calculation
        Eta_k = (alpha_k_n * eta_n) / (1 + ((alpha_k_n - 1) * eta_n))
        # Steam table calculation
        vk = Eta_k * fk
        lk = (1 - Eta_k) * fk
        Vtotal = np.sum(vk)
        Ltotal = np.sum(lk)
        xk = lk / Ltotal
        yk = vk / Vtotal
        print('\n The Flash calculation has convereged to the bubble point of majority component in liquid')
        print('\n The convergence temperature at BB point :- ', result)
        print('\n Stream properties:-')
        print('\n Vapour component flow rates :- ', vk)
        print('\n Liquid Component flow rates :- ', lk)
        print('\n Vapour component Mole Fraction :- ', xk)
        print('\n Liquid component Mole Fraction :- ', yk)
    elif x==False :
        P_spec = float(input("Enter the specified pressure "))
        temp_guess =float(input("Enter the Guess temperature "))
        result = fsolve(eq1_case_t_guess,temp_guess, args=(P_spec,A,B,C,fk,key,eta_n),xtol=1.7e-10)

        P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=result)
        Kn = P_vap / P_spec
        # calculation of the relative volatilities
        alpha_k_n = Kn / Kn[key - 1]
        # overhead split fraction calculation
        Eta_k = (alpha_k_n * eta_n) / (1 + ((alpha_k_n - 1) * eta_n))
        # Steam table calculation
        vk = Eta_k * fk
        lk = (1 - Eta_k) * fk
        Vtotal = np.sum(vk)
        Ltotal = np.sum(lk)
        xk = lk / Ltotal
        yk = vk / Vtotal
        print('\n The Flash calculation has convereged')
        print('\n The convergence temperature:- ', result)
        print ('\n Stream properties:-')
        print('\n Vapour component flow rates :- ',vk)
        print('\n Liquid Component flow rates :- ', lk)
        print('\n Vapour component Mole Fraction :- ', xk)
        print('\n Liquid component Mole Fraction :- ', yk)


def calculation(eta_n,P,key,fk,temp,A,B,C):
    P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=temp)
    Kn = P_vap / P
    # calculation of the relative volatilities
    alpha_k_n = Kn / Kn[key - 1]
    # overhead split fraction calculation
    Eta_k = (alpha_k_n * eta_n) / (1 + ((alpha_k_n - 1) * eta_n))
    # Steam table calculation
    vk = Eta_k * fk
    lk = (1 - Eta_k) * fk
    Vtotal = np.sum(vk)
    Ltotal = np.sum(lk)
    xk = lk / Ltotal
    temp_pcalc = (P_vap * xk)
    p_calc = np.sum(temp_pcalc)
    #print(Eta_k,eta_n)
    return (p_calc / P) - 1

def flash_case2(A,B,C,fk): #Flash P and T specified
    temp = float(input('Enter specified Temperature: '))
    P = float(input('Enter specified Pressure '))
    eta_n = float(input('Enter guess Eta for key component '))
    key = int(input("which component is the key"))
    P_vapcheck= vapour_pressure_temp(A=A, C=C, B=B, T=temp)
    P_bubble = np.sum(P_vapcheck*fk)
    P_dew = (np.sum(fk/P_vapcheck))**-1
    if P_bubble > P and P > P_dew: # Flash should be between Bubble and dew pressure
        print('\n Normal Calculation')
        result = fsolve(calculation, eta_n, args=(P,key,fk,temp,A,B,C), xtol=1e-10) #Fsolve method in scipy
        P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=temp)
        Kn = P_vap / P
        # calculation of the relative volatilities
        alpha_k_n = Kn / Kn[key - 1]
        # overhead split fraction calculation
        Eta_k = (alpha_k_n * result) / (1 + ((alpha_k_n - 1) * result))
        # Steam table calculation
        vk = Eta_k * fk
        lk = (1 - Eta_k) * fk
        Vtotal = np.sum(vk)
        Ltotal = np.sum(lk)
        xk = lk / Ltotal
        yk = vk / Vtotal
        print('\n The Flash calculation has convereged to the bubble point of majority component in liquid')
        print('\n The convergence En value is :- ', result)
        print('\n Stream properties:-')
        print('\n Vapour component flow rates :- ', vk)
        print('\n Liquid Component flow rates :- ', lk)
        print('\n Vapour component Mole Fraction :- ', xk)
        print('\n Liquid component Mole Fraction :- ', yk)
        print('\n Phi Value:- ', (Vtotal / np.sum(fk)))

        #print(Kn,alpha_k_n,Eta_k,P_vap  )

    else : print ('The system is not flashable')

def eq3_case_P_guess (P,temp,A,B,C,fk,phi,key): # if pressure is giving then this function is minimised
    P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=temp)
    # calculation of the partion coeeficients
    Kn = P_vap / P
    alpha_k_n = Kn / Kn[key - 1]
    theta = (Kn * phi) / (1 - phi)
    Eta_k = theta / (1 + theta)
    vk = Eta_k * fk
    lk = (1 - Eta_k) * fk
    Vtotal = np.sum(vk)
    Ltotal = np.sum(lk)
    xk = lk / Ltotal
    temp_pcalc = (P_vap * xk)
    p_calc = np.sum(temp_pcalc)
    # print(Eta_k,eta_n)
    return (p_calc / P) - 1

def eq3_case_T_guess (temp,P,A,B,C,fk,phi,key): # if temperature is specified then this function is minimised
    P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=temp)
    # calculation of the partion coeeficients
    Kn = P_vap / P
    alpha_k_n = Kn / Kn[key - 1]
    theta = (Kn*phi)/(1-phi)
    Eta_k= theta/(1+theta)
    vk = Eta_k * fk
    lk = (1 - Eta_k) * fk
    Vtotal = np.sum(vk)
    Ltotal = np.sum(lk)
    xk = lk / Ltotal
    temp_pcalc = (P_vap * xk)
    p_calc = np.sum(temp_pcalc)
    # print(Eta_k,eta_n)
    return (p_calc / P) - 1

def flash_case3(A,B,C,fk):# Case 3 where Phi and T/P are given
    key = int(input("which component is the key"))
    phi= float(input('Enter Specified Phi:-  '))
    ans=input('Are you going to specify pressure True/False?:- ')
    if ans.lower() =='true':
        x = True
    elif ans.lower() == 'false':
        x = False
    else:print ('Correct input was not given')
    if x == True:
        press = float(input('Enter the Specified Pressure:- '))
        temp = float(input('Enter Guess Temperature:-  '))
        result=fsolve(eq3_case_T_guess,temp,args=(press,A,B,C,fk,phi,key),xtol=1e-7)
        P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=result)
        # calculation of the partion coeeficients
        Kn = P_vap / press
        alpha_k_n = Kn / Kn[key - 1]
        theta = (Kn * phi) / (1 - phi)
        Eta_k = theta / (1 + theta)
        vk = Eta_k * fk
        lk = (1 - Eta_k) * fk
        Vtotal = np.sum(vk)
        Ltotal = np.sum(lk)
        xk = lk / Ltotal
        yk = vk / Vtotal

        return print('\n The Molar flow rate in Gas phases:- ',yk),print("\nThe molar flow in the liquid phase:- ",xk), print('\nThe converged temperature value:- ',result),
    elif x == False: #if pressure is not given the block is executed
        press = float(input('Enter the guess Pressure:- '))
        temp = float(input('Enter specified Temperature:-  '))
        result = fsolve(eq3_case_P_guess, press, args=(temp, A, B, C, fk, phi, key), xtol=1e-7)
        P_vap = vapour_pressure_temp(A=A, C=C, B=B, T=temp)
        # calculation of the partion coeeficients
        Kn = P_vap / result
        alpha_k_n = Kn / Kn[key - 1]
        theta = (Kn * phi) / (1 - phi)
        Eta_k = theta / (1 + theta)
        vk = Eta_k * fk
        lk = (1 - Eta_k) * fk
        Vtotal = np.sum(vk)
        Ltotal = np.sum(lk)
        xk = lk / Ltotal
        yk = vk / Vtotal
        return print('\nThe Molar flow rate in Gas phases:- ', yk), print("\nThe molar flow in the liquid phase:- ", xk), print(
            '\nThe converged pressure value:- ', result)



