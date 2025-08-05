from initial_conditions import InitialConditions, theta_prime
import numpy as np
from scipy.constants import G

# PLANET
# ------
rho_planet = 3.344 * 1e-3 / 1e-6  # Mittlere Dichte Mond (kg / m^3) (Wikipedia)
R_planet = 5700e3 # Radius Planet (m)
V_planet = 4./3. * np.pi * R_planet**3
M_planet = rho_planet * V_planet # Masse (kg)
T_planet =  24 * 60 * 60# Länge eines Tages (Dauer einer Umdrehung) (s)
omega_planet = 2.*np.pi/T_planet # Winkelgeschwindigkeit

# SATELIT
# -------
M_satelit = 235 # Masse (kg)
R_satelit = 2. * R_planet

# Koordinaten des Sateliten beim Zusammenstoss
phi_0 = 0.0
theta_0 = 0.0 # TODO: disen Parameter muss ich noch additiv anpassen
alpha = np.pi/4.  # Winkel der Anfangsgeschwindigkeit

# INITIAL CONIDTIONS
# ------------------

def get_initial_conditions(baffle: float) -> InitialConditions:

    v_satelit_vorher = np.sqrt( G*M_planet / R_satelit ) # (m/s)
    v_satelit_nachher = v_satelit_vorher * baffle
    theta_prime_0 = v_satelit_nachher / R_satelit

    return InitialConditions(
        theta_0=0.0,
        theta_prime_0=theta_prime_0,
        r_0=R_satelit,
        r_prime_0=0.0,
        m=M_satelit,
        M=M_planet)