from utils import compute_final_coordinates_newton
from configuration import get_initial_conditions_newton, Configuration
from utils import RootMemory, euclidean_to_spherical
from scipy.constants import G
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root


def main():
    # These are the actual coordinates of the cache
    # ---------------------------------------------
    cache_coordinate_N = 47.44670
    cache_coordinate_O = 8.12868
    # ---------------------------------------------

    initial_conditions_newton = get_initial_conditions_newton()

    root_memory = RootMemory()

    def gravity_rhs(
            t: float,
            u: tuple[float, float, float, float, float, float]) -> tuple[float, float, float, float, float, float]:
        """
        returns the rhs of the gravitational ODE, assuming that
        u in R^6 is the euclidean coordinate and velocity, i.e.
        u = [x, y, z, x_dot, y_dot, z_dot]
        """
        k = G * Configuration.M_planet
        x, y, z, x_dot, y_dot, z_dot = u

        r = np.sqrt(x**2 + y**2 + z**2)

        return [x_dot, y_dot, z_dot, -k/r**3 * x, -k/r**3 * y, -k/r**3 * z]

    def F(t_end):
        result = solve_ivp(
            fun=gravity_rhs,
            t_span=(0.0, t_end),
            y0=initial_conditions_newton.Y_0, rtol=1e-8)
        u = result.y
        t = result.t

        r_end = np.sqrt(u[0, -1]**2 + u[1, -1]**2 + u[2, -1]**2)

        root_memory.add_to_memory(
            t=t[-1], u=u[:, -1])

        return r_end - Configuration.R_planet

    _ = root(fun=F, x0=10.0, options={'xtol': 1e-8}).x

    t_hit = root_memory.ts[-1]

    u_hit = root_memory.us[-1]
    x_hit, y_hit, z_hit = u_hit[0:3]

    _, theta_hit, phi_hit = euclidean_to_spherical(
        x=x_hit, y=y_hit, z=z_hit)

    # looking at the definition of spherical coordinates,
    # we need to transform the theta coordinate to receive
    # the angle measured from equator instead of the angle
    # measured from north pole
    theta_hit_tilde = np.pi/2. - theta_hit

    # planet is rotating during the fall
    phi_hit -= t_hit * Configuration.omega_planet

    O = (phi_hit + Configuration.phi_0) / (2.*np.pi) * 360.
    N = theta_hit_tilde / (2.*np.pi) * 360.

    print(f'GeoCache:   (N, O) \t= ({cache_coordinate_N:.5f}, \t\t{cache_coordinate_O:.5f})')
    print(f'Experiment: (N, O) \t= ({N}, \t{O})')


if __name__ == '__main__':
    main()
