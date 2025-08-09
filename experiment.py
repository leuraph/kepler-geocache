
from configuration import get_initial_conditions
from utils import compute_theta_hit_ode, compute_theta_hit_integration, compute_theta_hit_closed_form
from utils import compute_dt_from_ode_via_theta_hit, compute_dt_from_ode_via_radius


def main() -> None:
    baffle = 0.565

    # EXPERIMENT
    # ----------
    initial_conditions = get_initial_conditions(baffle=baffle)

    # THETA HIT
    # ---------
    print(f'Computing theta_hit using different approaches...')

    theta_hit_closed_form = compute_theta_hit_closed_form(
        initial_conditions=initial_conditions)
    theta_hit_integration = compute_theta_hit_integration(initial_conditions=initial_conditions)[0]

    # USING ODE
    t_hit_from_radius = compute_dt_from_ode_via_radius(initial_conditions=initial_conditions)
    t_hit_from_theta_hit_closed = compute_dt_from_ode_via_theta_hit(initial_conditions=initial_conditions, theta_hit=theta_hit_closed_form)
    t_hit_from_theta_hit_integation = compute_dt_from_ode_via_theta_hit(initial_conditions=initial_conditions, theta_hit=theta_hit_integration)

    print(f'dt 1 = {t_hit_from_radius}')
    print(f'dt 1 = {t_hit_from_theta_hit_closed}')
    print(f'dt 1 = {t_hit_from_theta_hit_integation}')
    
    theta_hit_ode_1 = compute_theta_hit_ode(initial_conditions=initial_conditions, t_hit=t_hit_from_radius)
    theta_hit_ode_2 = compute_theta_hit_ode(initial_conditions=initial_conditions, t_hit=t_hit_from_theta_hit_closed)
    theta_hit_ode_3 = compute_theta_hit_ode(initial_conditions=initial_conditions, t_hit=t_hit_from_theta_hit_integation)

    print(f'theta_hit closed form = {theta_hit_closed_form}')
    print(f'theta_hit integration = {theta_hit_integration}')
    print(f'theta_hit ode 1 = {theta_hit_ode_1}')
    print(f'theta_hit ode 2 = {theta_hit_ode_2}')
    print(f'theta_hit ode 3 = {theta_hit_ode_3}')


if __name__ == '__main__':
    main()
