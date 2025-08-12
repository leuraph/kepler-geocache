from mimetypes import init
import numpy as np
from configuration import R_planet, R_satelit, get_initial_conditions, omega_planet
import initial_conditions
from scipy.integrate import solve_ivp
from scipy.optimize import root
from utils import euclidean_to_spherical, gravity_rhs


class RootMemory:
    ts: list[float]
    us: list[tuple[float, float, float, float, float, float]]

    def __init__(self):
        self.ts = []
        self.us = []

    def add_to_memory(self, t: float, u: tuple[float, float, float, float, float, float]):
        self.ts.append(t)
        self.us.append(u)


def main():
    baffle = 0.6

    # provide the initial conditions
    # u_0 = [x_0, y_0, z_0, x_dot_0, y_dot_0, z_dot_0]
    initial_conditions = get_initial_conditions(baffle=baffle)
    # Coordinates
    x_0 = R_satelit
    y_0 = 0.0
    z_0 = 0.0
    # velocities
    x_dot_0 = 0
    y_dot_0 = initial_conditions.v_0 / np.sqrt(2.0)
    z_dot_0 = initial_conditions.v_0 / np.sqrt(2.0)

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
    print(t_hit)

    u_hit = root_memory.us[-1]
    print(u_hit)
    x_hit, y_hit, z_hit = u_hit[0:3]

    r_hit, theta_hit, phi_hit = euclidean_to_spherical(
        x=x_hit, y=y_hit, z=z_hit)
    
    theta_hit_tilde = np.pi/2. - theta_hit
    
    d_angle = np.arccos(np.cos(phi_hit)*np.cos(theta_hit_tilde))
    beta = np.arcsin( np.sin(theta_hit_tilde) / np.sin(d_angle))

    print(f'd_angle = {d_angle}')
    print(f'beta = {beta}, pi/4 = {np.pi/4.}')

    # planet is rotating during the fall
    phi_hit -= t_hit * omega_planet

    print(f'phi_hit = {phi_hit}')

    O = phi_hit / (2.*np.pi) * 360.
    N = theta_hit_tilde / (2.*np.pi) * 360.

    print(f'(N, O) = ({N}, {O})')

    

if __name__ == '__main__':
    main()
