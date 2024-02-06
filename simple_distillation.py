import math
import numpy as np



def antoine_eqn(T, A, B, C):
    return 10**(A - B / (T + C))
def antonie_eq_e (T, A, B, C):
    return math.e**(A - B / (T + C))
def temp_vapour_pressure (A,B,C,P):
    return (B/(A-math.log(P)))-C

def distillationtower(A,B,C,fk,T,P,condenser_check=False):
    heavy_key_comp = int (input ('Enter the which is a heavy key:- '))
    light_key_comp = int (input ('\nEnter which is the light key:- '))
    heavy_key_split = float (input ('\nEnter the Heavy Key split'))
    light_key_split = float (input ('\nEnter the light key split'))
    ans= input ('\n Is there a distributed component:- ')

    P_vap= antonie_eq_e(T,A,B,C) #Vapour Pressure Array
    # partition coefficent
    Kn = P_vap / P
    # calculation of the relative volatilities
    alpha_k_n = Kn / Kn[heavy_key_comp - 1] # Relative volatilities calculated on the basis of
    #Fenske Eq stage calculation minimum trays
    N= (math.log((light_key_split*(1+heavy_key_split))/(heavy_key_split*(1-light_key_split))))/math.log(alpha_k_n[light_key_comp-1]/alpha_k_n[heavy_key_comp-1])
    N_m = math.ceil(N) # Number of Stages Final
    #distributed component numpy array

    if heavy_key_comp-light_key_comp == 1:
        print('\nNo distributed components !')
    else:
        for i in range (light_key_comp+1,heavy_key_comp):
            eta_k= (((alpha_k_n[i-1])**N_m)*heavy_key_split)/(1+(((alpha_k_n[i-1])**N_m)-1)*heavy_key_split)
            if eta_k.size == 0:
                print('\nError in calculation.')
            else: print('\n All good. ')
    print ('\nThe number of stages is :- ', N_m)
    eta_k_list = [light_key_split,heavy_key_split]
    insert_index = 1
    eta_k_final = eta_k_list[:insert_index] + [eta_k.tolist()] + eta_k_list[insert_index:]
    print ('The Split fraction based on tops is  :- ',eta_k_final)
    distilate = eta_k_final * fk

    bottoms = (np.full(np.size(eta_k_final),1)-eta_k_final) * fk
    total_bottoms = np.sum(bottoms)
    total_distillate = np.sum(distilate)
    mole_frac_tops = distilate /total_distillate
    mole_frac_bottoms = bottoms /total_bottoms
    print ('\n The distillate flows are :- ', distilate)
    print('\n The distillate mole fractions  are :- ', mole_frac_tops)
    print('\n The bottoms flows are :- ', bottoms)
    print('\n The bottoms mole fractions are :- ', distilate)
    if condenser_check:
        P_vap_bubble = antonie_eq_e(T, A, B, C)  # Vapour Pressure Array
        # partition coefficent
        Kn_bubble = P_vap / P
        # calculation of the relative volatilities
        alpha_k_n = Kn / Kn[light_key_comp-1]
        alpha_bar_bubble = np.sum(mole_frac_tops*alpha_k_n)
        P_o= P/alpha_bar_bubble
        print("The Bubble point tempearature of the most volatile content in the Distilate:- ",temp_vapour_pressure(A=A[light_key_comp-1],B=B[light_key_comp-1],C=C[light_key_comp-1],P=P_o))






'''
A=[15.9008,16.0137,16.1156]
Ai=np.array(A)
B =[2788.51,3096.52,3395.57]
Bi=np.array(B)
C =[-52.34,-53.67,-59.44]
Ci=np.array(C)
fi=np.array([20,30,50])

distillationtower(A,B,C,fi,T=np.array([386]),P=1,condenser_check=False)

'''

