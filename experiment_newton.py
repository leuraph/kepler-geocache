from utils import compute_final_coordinates_newton


def main():
    # These are the actual coordinates of the cache
    # ---------------------------------------------
    cache_coordinate_N = 47.44670
    cache_coordinate_O = 8.12868
    # ---------------------------------------------

    baffle = 0.53283236
    O_initial = -5.287518

    N, O = compute_final_coordinates_newton(
        baffle=baffle, O_initial=O_initial)

    print(f'(N, O) = ({N}, {O})')
    print(f'(N_cache, O_cache) = ({cache_coordinate_N}, {cache_coordinate_O})')


if __name__ == '__main__':
    main()
