import math

# Constants
k_B = 1.38e-23  # Boltzmann constant in J/K
h = 6.626e-34   # Planck's constant in J*s
T = 298         # Temperature in Kelvin
sigma_CH3=3
sigma_C2H6 =18
sigma_coll = 0.46e-18  # Collision cross-section in m^2
c= 3e10
# CH3 data
m_CH3 = 2.49081e-26  # Mass of CH3 in kg
B_CH3_cm = 14.0     # Rotational constant for CH3 in m^-1

# C2H6# data
m_C2H6 = 4.98151e-26  # Mass of C2H6# in kg
B_C2H6_cm = 2.80      # Rotational constant for C2H6# in cm^-1

# Convert rotational constants to m^-1
A_CH3_m = 9.32 * 100
A_C2H6_m = 4.66 * 100
B_CH3_m = B_CH3_cm * 100
B_C2H6_m = B_C2H6_cm * 100
C_CH3_m = B_CH3_cm * 100
C_C2H6_m = B_C2H6_cm * 100
# Translational contributions
q_trans_CH3 = (2 * 3.14 * m_CH3 * k_B * T / h**2)**1.5
q_trans_C2H6 = (2 * 3.14 * m_C2H6 * k_B * T / h**2)**1.5
print('CH3 :- Translation :',q_trans_CH3)
print('C2H6#:- Translation :',q_trans_C2H6)
# Rotational contributions
q_rot_CH3 = (1/sigma_CH3)*(((k_B*T)/(h*c))**1.5)*((3.14/(A_CH3_m *B_CH3_m*C_CH3_m))**0.5)
q_rot_C2H6 =(1/sigma_C2H6)*(((k_B*T)/(h*c))**1.5)*((3.14/(A_C2H6_m *B_C2H6_m*C_C2H6_m))**0.5)
# non-linear moleculecules
print('CH3 :- Rotational :',q_trans_CH3)
print('C2H6#:- Rotational :',q_trans_C2H6)
# Vibrational contributions
q_vib_CH3 = 1
q_vib_C2H6 = 1

# Electronic contributions
q_elec_CH3 = 1
q_elec_C2H6 = 1

# Molecular partition functions
Q_CH3 = q_trans_CH3 * q_rot_CH3 * q_vib_CH3 * q_elec_CH3
Q_C2H6 = q_trans_C2H6 * q_rot_C2H6 * q_vib_C2H6 * q_elec_C2H6
print('The total molecular partition Function Ch3:- ',Q_CH3)
print('The total molecular partition Function C2H6:- ',Q_C2H6)
# Pre-exponential factors from Transition State Theory
A_TST_CH3 = (k_B * T / h) / Q_CH3
A_TST_C2H6 = (k_B * T / h) / Q_C2H6

# Pre-exponential factor from Collision Theory
mu = m_CH3 * m_C2H6 / (m_CH3 + m_C2H6)  # Reduced mass
A_coll = sigma_coll * math.sqrt(k_B * T / (2 * 3.14 * mu))

# Steric Factor
steric_factor = A_TST_CH3 / A_coll

# Display results
print(f"Pre-exponential factor (TST) for CH3: {A_TST_CH3:.3e} s^-1")
print(f"Pre-exponential factor (TST) for C2H6#: {A_TST_C2H6:.3e} s^-1")
print(f"Pre-exponential factor (Collision Theory): {A_coll:.3e} m^3/(mol*s)")
print(f"Steric Factor: {steric_factor:.3e}")
