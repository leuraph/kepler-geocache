from scipy.optimize import root
from configuration import get_initial_conditions
from scipy.integrate import solve_ivp
from scipy.constants import G
import matplotlib.pyplot as plt
from configuration import R_planet, R_satelit, alpha, phi_0, theta_0, omega_planet
import numpy as np
from scipy.integrate import quad


def main() -> None:
    baffle = 0.565

    # EXPERIMENT
    # ----------
    initial_conditions = get_initial_conditions(baffle=baffle)

    def kepler_rhs(t: float, y: tuple[float, float]) -> tuple[float, float]:
        r, r_dot, theta = y
        theta_dot = initial_conditions.get_conserved_angular_momentum() / (initial_conditions.m * r**2)
        r_ddot = r * theta_dot**2 - G * initial_conditions.M / r**2
        return [r_dot, r_ddot, theta_dot]

    y0 = [initial_conditions.r_0, initial_conditions.r_prime_0, initial_conditions.theta_0]
    t_guess = 1.0 # (s)

    # ZEITPUNKT DES EINSCHLAGS
    # ------------------------
    def F(t_end):

        result = solve_ivp(fun=kepler_rhs, t_span=(0., t_end), y0=y0, rtol=1e-8)
        t = result.t
        y = result.y
        r_end = y[0, -1]

        return r_end - R_planet

    t_hit = root(fun=F, x0=t_guess, options={'xtol': 0.000001}).x
    print(f'Satelit schlägt ein nach {t_hit / 3600} h')

    # WINKEL BEIM EINSCHLAG
    # ---------------------
    t_end = t_hit
    result = solve_ivp(fun=kepler_rhs, t_span=(0., t_hit), y0=y0, rtol=1e-8)

    theta_hit = result.y[2, -1]
    print(f'Überflogener Winkel bis zum Einschlag {theta_hit} (Bogenmass)')

    # Umrechnung in Koordinaten
    def get_final_coordinates(theta_hit) -> tuple[float, float]:
        d_theta = np.arcsin( np.sin(theta_hit)  * np.sin(alpha))
        d_phi = np.arccos( np.cos(theta_hit) / np.cos(d_theta) )

        phi = phi_0 + d_phi # Kollision über Äquator
        theta = theta_0 + d_theta - t_hit * omega_planet

        N = 360. * phi / (2.* np.pi)  # Breitengrad
        O = 360. * theta / (2. * np.pi)  # Längengrad
        return N, O

    N, O = get_final_coordinates(theta_hit=theta_hit)
    print(f'Einschlag bei (N, O) = ({N}, {O})')

    print(f'Sanity check mit analytischer Lösungsmethode + Quadratur:')
    def integrand(r):
        l = initial_conditions.get_conserved_angular_momentum()
        E = initial_conditions.get_conserved_total_energy()
        V_eff = l**2 / (2 * initial_conditions.m * r**2) - G*initial_conditions.M*initial_conditions.m/r
        radicand = 2 * initial_conditions.m * (E - V_eff)
        return l / (r**2 * np.sqrt(radicand))
    
    theta_hit = quad(func=integrand, a=R_planet, b=R_satelit)
    print(f'Überflogener Winkel bis zum Einschlag {theta_hit} (Bogenmass)')
    N, O = get_final_coordinates(theta_hit=theta_hit)
    print(f'Einschlag bei (N, O) = ({N}, {O})')

if __name__ == '__main__':
    main()
