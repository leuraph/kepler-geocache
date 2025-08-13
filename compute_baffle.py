from utils import compute_final_coordinates_newton
from scipy.optimize import root


def main() -> None:

    cache_coordinate_N = 47.44670
    cache_coordinate_O = 8.12868

    baffle_0 = 0.53283236

    def F_baffle(baffle: float):
        print(float(baffle))
        N, O = compute_final_coordinates_newton(
            baffle=float(baffle), O_initial=0.0)
        print(f'N = {N}')
        return N - cache_coordinate_N
    
    baffle = root(fun=F_baffle, x0=baffle_0, options={'xtol': 1e-10}).x

    print(f'baffle = {baffle}')


if __name__ == '__main__':
    main()
