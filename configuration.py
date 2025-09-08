from dataclasses import dataclass
from initial_conditions import InitialConditions, angular_momentum
import numpy as np
from scipy.constants import G
from physics import get_total_energy, get_angular_momentum

@dataclass
class Configuration:
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
    gamma = 0.53283236

    # Koordinaten des Sateliten beim Zusammenstoss
    phi_0 = -5.287518 * (2.*np.pi/360.)  # (West / Ost)
    theta_0 = 0.0  # (Süd / Nord)
    alpha = np.pi/8.  # Winkel der Anfangsgeschwindigkeit


@dataclass
class InitialConditionsNewton:
    initial_position: tuple[float, float, float]  # euclidean initial coordinates
    initial_velocity: tuple[float, float, float]  # euclidean initial velocity
    Y_0: tuple[float, float, float, float, float, float]



def get_initial_conditions_newton() -> InitialConditionsNewton:
    # position
    x_0 = R_satelit
    y_0 = 0.0
    z_0 = 0.0
    u_0 = [x_0, y_0, z_0]

    v_0 = np.sqrt( G*Configuration.M_planet / Configuration.R_satelit ) # (m/s)
    v_0_tilde = v_0 * Configuration.gamma  # velocity after collision

    # velocity components
    x_dot_0 = 0.
    y_dot_0 = np.sin(Configuration.alpha) * v_0_tilde
    z_dot_0 = np.cos(Configuration.alpha) * v_0_tilde
    u_dot_0 = [x_dot_0, y_dot_0, z_dot_0]

    # combined 6-dimensional components
    Y_0 = [x_0, y_0, z_0, x_dot_0, y_dot_0, z_dot_0]

    return InitialConditionsNewton(
        initial_position=u_0,
        initial_velocity=u_dot_0,
        Y_0=Y_0)


@dataclass
class InitialConditionsKepler:
    # initial polar coordinates
    theta_0: float
    r_0: float

    # initial velocity in polar coordinates
    theta_0_dot: float
    r_0_dot: float

    # conserved quantities
    angular_momentum: float
    total_energy: float


def get_initial_conditions_kepler() -> InitialConditionsKepler:
    theta_0 = 0.0
    r_0 = Configuration.R_satelit

    v_0 = np.sqrt( G*Configuration.M_planet / Configuration.R_satelit ) # (m/s)
    v_0_tilde = v_0 * Configuration.gamma  # velocity after collision

    theta_0_dot = v_0_tilde / Configuration.R_satelit
    r_0_dot = 0.0

    angular_momentum = get_angular_momentum(
        m=Configuration.M_satelit,
        r=Configuration.R_satelit,
        theta_dot=theta_0_dot)

    total_energy = get_total_energy(
        m=Configuration.M_satelit,
        M=Configuration.M_planet,
        r=Configuration.R_satelit,
        r_dot=r_0_dot,
        l=angular_momentum)

    return InitialConditionsKepler(
        theta_0=theta_0,
        r_0=r_0,
        theta_0_dot=theta_0_dot,
        r_0_dot=r_0_dot,
        angular_momentum=angular_momentum,
        total_energy=total_energy)


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
alpha = np.pi/8.  # Winkel der Anfangsgeschwindigkeit

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
        M=M_planet,
        alpha=alpha)
