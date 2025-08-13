from utils import compute_final_coordinates_newton
from scipy.optimize import root


def main() -> None:

    cache_coordinate_N = 47.44670
    cache_coordinate_O = 8.12868

    baffle = 0.53283236

    def F_east(O_initial: float):
        print(float(O_initial))
        N, O = compute_final_coordinates_newton(
            baffle=baffle, O_initial=float(O_initial))
        print(f'O = {O}')
        return O - cache_coordinate_O
    
    O_initial = root(fun=F_east, x0=-5.28751838, options={'xtol': 1e-10}).x

    print(f'O_initial = {O_initial}')


if __name__ == '__main__':
    main()
