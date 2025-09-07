from initial_conditions import InitialConditions, effective_potential
from configuration import R_planet, R_satelit, get_initial_conditions, omega_planet
from scipy.optimize import root
from scipy.integrate import quad
import numpy as np
from scipy.constants import G
from scipy.integrate import solve_ivp



def theta_prime(m: float, l: float, r: float) -> float:
    return l / (m * r**2)

def r_prime(m: float, M: float, r: float, l: float, E: float) -> float:
    return np.sqrt(
        2. / m *(
            E - effective_potential(m=m, M=M, r=r, l=l)
        )
    )


def get_final_coordinates(
        theta_initial, phi_initial,
        theta_hit: float, t_hit,
        omega_planet,
        initial_conditions: InitialConditions) -> tuple[float, float]:
    """trigonmetry on a sphere"""
    beta = np.pi/2. - initial_conditions.alpha
    d_theta = np.arcsin( np.sin(theta_hit) * np.sin(beta) )
    d_phi = np.arccos(np.cos(theta_hit) / np.cos(d_theta))

    phi = phi_initial + d_phi - t_hit * omega_planet
    theta = theta_initial + d_theta

    N = theta * 360. / (2.*np.pi)  # Breitengrad
    O = phi * 360. / (2.*np.pi)  # Längengrad
    return N, O


def kepler_rhs(
        initial_conditions: InitialConditions,
        t: float,
        y: tuple[float, float, float]) -> tuple[float, float, float]:
    r, r_dot, theta = y
    theta_dot = initial_conditions.conserved_angular_momentum / (initial_conditions.m * r**2)
    r_ddot = r * theta_dot**2 - G * initial_conditions.M / r**2
    return [r_dot, r_ddot, theta_dot]


def gravity_rhs(
        initial_conditions: InitialConditions,
        t: float,
        u) -> tuple[float, float, float, float, float, float]:
    """
    returns the rhs of the gravitational ODE, assuming that
    u in R^6 is the euclidean coordinate and velocity, i.e.
    u = [x, y, z, x_dot, y_dot, z_dot]
    """
    k = G * initial_conditions.M
    x, y, z, x_dot, y_dot, z_dot = u

    r = np.sqrt(x**2 + y**2 + z**2)

    return [x_dot, y_dot, z_dot, -k/r**3 * x, -k/r**3 * y, -k/r**3 * z]


def spherical_to_euclidean(r: float, theta: float, phi: float) -> tuple[float, float, float]:
    """
    Given spherical coordinates (r, theta, phi), returns Euclidean (x, y, z)
    theta: polar angle from +z axis
    phi: azimuthal angle from +x axis in xy-plane
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def euclidean_to_spherical(x: float, y: float, z: float) -> tuple[float, float, float]:
    """
    Given Euclidean coordinates (x, y, z), returns spherical (r, theta, phi)
    r: radial distance
    theta: polar angle from +z axis
    phi: azimuthal angle from +x axis in xy-plane
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    if r == 0:
        return 0.0, 0.0, 0.0  # Convention for origin

    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return r, theta, phi


def eccentricity(
        total_energy: float,
        m: float,
        M: float,
        angular_momentum: float) -> float:
    """source: https://de.wikipedia.org/wiki/Zweik%C3%B6rperproblem"""
    h = angular_momentum / m
    e = np.sqrt(1. + (2. * total_energy * h**2)/((G*M)**2 * m))
    return e


def radius(
        theta: float,
        theta_0: float,
        initial_conditions: InitialConditions) -> float:
    """https://de.wikipedia.org/wiki/Zweik%C3%B6rperproblem"""
    m = initial_conditions.m
    M = initial_conditions.M
    angular_momentum = initial_conditions.conserved_angular_momentum
    total_energy = initial_conditions.conserved_total_energy

    h = angular_momentum / m
    e = eccentricity(
        total_energy=total_energy,
        m=m, M=M, angular_momentum=angular_momentum)
    numerator = h**2 / (G*M)
    denominator = 1. + e * np.cos(theta - theta_0)
    return numerator / denominator


def compute_theta_hit_closed_form(initial_conditions: InitialConditions) -> float:
    """source: https://de.wikipedia.org/wiki/Zweik%C3%B6rperproblem"""
    def F(theta):
        final_radius = radius(
            theta=theta,
            theta_0=np.pi,
            initial_conditions=initial_conditions)
        return final_radius - R_planet
    theta_hit = root(fun=F, x0=1.0, options={'xtol': 1e-8}).x
    return theta_hit


def compute_theta_hit_integration(
        initial_conditions: InitialConditions) -> float:
    def integrand(r):
        l = initial_conditions.get_conserved_angular_momentum()
        E = initial_conditions.get_conserved_total_energy()
        V_eff = l**2 / (2 * initial_conditions.m * r**2) - G*initial_conditions.M*initial_conditions.m/r
        radicand = 2 * initial_conditions.m * (E - V_eff)
        return l / (r**2 * np.sqrt(radicand))
    
    theta_hit = quad(func=integrand, a=R_planet, b=R_satelit)
    return theta_hit


def compute_theta_hit_ode(initial_conditions: InitialConditions, t_hit: float) -> float:

    y0 = [initial_conditions.r_0, initial_conditions.r_prime_0, initial_conditions.theta_0]

    def rhs(t, y):
        return kepler_rhs(
            initial_conditions=initial_conditions, t=t, y=y)

    result = solve_ivp(fun=rhs, t_span=(0., t_hit), y0=y0, rtol=1e-8)
    y = result.y
    theta_hit = y[2, -1]
    return theta_hit


def compute_dt_from_ode_via_theta_hit(
        initial_conditions: InitialConditions,
        theta_hit) -> float:

    y0 = [initial_conditions.r_0, initial_conditions.r_prime_0, initial_conditions.theta_0]
    t_guess = 1.0 # (s)

    def rhs(t, y):
        return kepler_rhs(
            initial_conditions=initial_conditions, t=t, y=y)

    # ZEITPUNKT DES EINSCHLAGS
    # ------------------------
    def F(dt):

        result = solve_ivp(fun=rhs, t_span=(0., dt), y0=y0, rtol=1e-8)
        t = result.t
        y = result.y
        theta_end = y[2, -1]

        return theta_end - theta_hit

    t_hit = root(fun=F, x0=t_guess, options={'xtol': 1e-8}).x
    return t_hit


def compute_dt_from_ode_via_radius(
        initial_conditions: InitialConditions) -> float:

    y0 = [
        initial_conditions.r_0,
        initial_conditions.r_prime_0,
        initial_conditions.theta_0]
    t_guess = 1.0 # (s)

    def rhs(t, y):
        return kepler_rhs(
            initial_conditions=initial_conditions, t=t, y=y)

    # ZEITPUNKT DES EINSCHLAGS
    # ------------------------
    def F(dt):

        result = solve_ivp(fun=rhs, t_span=(0., dt), y0=y0, rtol=1e-8)
        t = result.t
        y = result.y
        r_end = y[0, -1]

        return r_end - R_planet

    t_hit = root(fun=F, x0=t_guess, options={'xtol': 1e-8}).x
    return t_hit


def compute_dt_from_ode_via_theta_hit(
        initial_conditions: InitialConditions,
        theta_hit: float) -> float:

    y0 = [
        initial_conditions.r_0,
        initial_conditions.r_prime_0,
        initial_conditions.theta_0]
    t_guess = 1.0 # (s)

    def rhs(t, y):
        return kepler_rhs(
            initial_conditions=initial_conditions, t=t, y=y)

    # ZEITPUNKT DES EINSCHLAGS
    # ------------------------
    def F(dt):

        result = solve_ivp(fun=rhs, t_span=(0., dt), y0=y0, rtol=1e-8)
        t = result.t
        y = result.y
        theta_end = y[2, -1]

        return theta_end - theta_hit

    t_hit = root(fun=F, x0=t_guess, options={'xtol': 1e-8}).x
    return t_hit

class RootMemory:
    ts: list[float]
    us: list[tuple[float, float, float, float, float, float]]

    def __init__(self):
        self.ts = []
        self.us = []

    def add_to_memory(self, t: float, u: tuple[float, float, float, float, float, float]):
        self.ts.append(t)
        self.us.append(u)


def compute_final_coordinates_newton(
        baffle: float,
        O_initial: float) -> tuple[float, float]:

    # provide the initial conditions
    # u_0 = [x_0, y_0, z_0, x_dot_0, y_dot_0, z_dot_0]
    initial_conditions = get_initial_conditions(baffle=baffle)
    # Coordinates
    x_0 = R_satelit
    y_0 = 0.0
    z_0 = 0.0
    # velocities
    x_dot_0 = 0
    y_dot_0 = np.sin(initial_conditions.alpha) * initial_conditions.v_0
    z_dot_0 = np.cos(initial_conditions.alpha) * initial_conditions.v_0

    u_0 = [x_0, y_0, z_0, x_dot_0, y_dot_0, z_dot_0]

    root_memory = RootMemory()

    def gravity_rhs_fun(t, y):
        return gravity_rhs(
            initial_conditions=initial_conditions, t=t, u=y)

    def F(t_end):
        result = solve_ivp(
            fun=gravity_rhs_fun,
            t_span=(0.0, t_end),
            y0=u_0, rtol=1e-8)
        u = result.y
        t = result.t

        r_end = np.sqrt(u[0, -1]**2 + u[1, -1]**2 + u[2, -1]**2)

        root_memory.add_to_memory(
            t=t[-1], u=u[:, -1])

        return r_end - R_planet

    _ = root(fun=F, x0=10.0, options={'xtol': 1e-8}).x

    t_hit = root_memory.ts[-1]

    u_hit = root_memory.us[-1]
    x_hit, y_hit, z_hit = u_hit[0:3]

    _, theta_hit, phi_hit = euclidean_to_spherical(
        x=x_hit, y=y_hit, z=z_hit)
    
    theta_hit_tilde = np.pi/2. - theta_hit
    
    d_angle = np.arccos(np.cos(phi_hit)*np.cos(theta_hit_tilde))
    beta = np.arcsin( np.sin(theta_hit_tilde) / np.sin(d_angle))

    print(f'd_angle = {d_angle}')
    print(f'beta_comp = {beta}, beta_initial = {np.pi/2. - initial_conditions.alpha}')

    # planet is rotating during the fall
    phi_hit -= t_hit * omega_planet

    O = phi_hit / (2.*np.pi) * 360. + O_initial
    N = theta_hit_tilde / (2.*np.pi) * 360.

    return N, O

def get_angular_momentum(
        m: float,
        r: float,
        theta_dot: float) -> float:
    return m * r**2 * theta_dot

def get_potential(m: float, M: float, r:float) -> float:
    return - G * M * m / r

def get_total_energy(
        m: float,
        M: float,
        r:float,
        r_dot: float,
        l:float) -> float:
    return (
        0.5 * m * r_dot ** 2
        + get_effective_potential(m=m, M=M, r=r, l=l)
        )

def get_effective_potential(
        m: float,
        M: float,
        r: float,
        l: float) -> float:
    return get_potential(m=m, M=M, r=r) + l**2 / (2. * m * r**2)