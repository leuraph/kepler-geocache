from scipy.constants import G

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