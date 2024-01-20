import numpy as np
import math
from scipy.optimize import fsolve

# creating a function for Vapour pressure calculation from temperature
def vapour_pressure_temp (A,B,C,T):
    return math.e**(A - (B / (T + C))) #Antonie eq in e form
# creating a function to return temperature from Vapor Pressure
def temp_vapour_pressure (A,B,C,P):
    return (B/(A-math.log(P,math.e)))-C
# making a function for relative volitality calculations
#making function for Case one of Flash calculation
def flash_case1(A,B,C,fk):
    eta_n = float(input('Enter the overhead split fraction of the key component:  '))
    x= input("Specify Temperature yes/no?  ")
    key= int(input("which component is the key"))
    def eq1_case_P_guess(temp, P):
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
        yk = vk / Vtotal
        # finding the highest liquid phase component
        index = np.argmax(lk)
        # Alpha bar calculation
        alpha_bar = np.sum(alpha_k_n * xk)
        P_calc = (alpha_bar/alpha_k_n[index])*P_vap[index]
        return P_calc - P_vap[index]

    def eq1_case_t_guess(temp, P):
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
        yk = vk/Vtotal
        # finding the highest liquid phase component
        index = np.argmax(lk)
        # Alpha bar calculation
        alpha_bar = np.sum(alpha_k_n * xk)
        P_calc = (alpha_k_n[index]/alpha_bar)*P
        return P_calc - P_vap[index]
    if x.lower() == 'yes':
        Temp_spec= float(input("Enter the specified temperature in Kelvin: "))
        P_guess1 = float(input("Enter the guess Pressure in mmHg"))
        result= fsolve(eq1_case_P_guess,P_guess1,args=Temp_spec)
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
        print('\n The comvergence temperature at BB point :- ', result)
        print('\n Stream properties:-')
        print('\n Vapour component flow rates :- ', vk)
        print('\n Liquid Component flow rates :- ', lk)
        print('\n Vapour component Mole Fraction :- ', xk)
        print('\n Liquid component Mole Fraction :- ', yk)
    elif x.lower() == 'no':
        P_spec = float(input("Enter the specified pressure in mmHg"))
        temp_guess =float(input("Enter the Guess temperature in Kelvin "))
        result = fsolve(eq1_case_t_guess, temp_guess, args=(P_spec))
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
        print('\n The Flash calculation has convereged to the bubble point of majority component in liquid')
        print('\n The comvergence temperature at BB point :- ', result)
        print ('\n Stream properties:-')
        print('\n Vapour component flow rates :- ',vk)
        print('\n Liquid Component flow rates :- ', lk)
        print('\n Vapour component Mole Fraction :- ', xk)
        print('\n Liquid component Mole Fraction :- ', yk)


A =[15.9008,16.0137,16.1156]
Ai=np.array(A)
B =[2788.51,3096.52,3395.57]
Bi=np.array(B)
C =[-52.34,-53.67,-59.44]
Ci=np.array(C)
fi=np.array([30,50,40])
flash_case1(Ai,Bi,Ci,fi)