import numpy as np
from scipy.constants import G


class InitialConditions:
    theta_0: float
    theta_prime_0: float
    r_0: float
    r_prime_0: float
    m: float
    M: float
    k: float

    def __init__(
            self,
            theta_0: float,
            theta_prime_0: float,
            r_0: float,
            r_prime_0: float,
            m: float,
            M: float):
        self.theta_0 = theta_0
        self.theta_prime_0 = theta_prime_0
        self.r_0 = r_0
        self.r_prime_0 = r_prime_0
        self.m = m
        self.M = M
        self.k = G * self.M * self.m
    
    def get_conserved_angular_momentum(self) -> float:
        return angular_momentum(
            m=self.m,
            r=self.r_0,
            theta_prime=self.theta_prime_0)
    
    def get_conserved_total_energy(self) -> float:
        return total_energy(
            m=self.m,
            M=self.M,
            r=self.r_0,
            r_prime=self.r_prime_0,
            l=self.get_conserved_angular_momentum())


def angular_momentum(
        m: float,
        r: float,
        theta_prime: float) -> float:
    return m * r**2 * theta_prime

def potential(m: float, M: float, r:float) -> float:
    return - G * M * m / r

def total_energy(
        m: float,
        M: float,
        r:float,
        r_prime: float,
        l:float) -> float:
    return (
        0.5 * m * r_prime ** 2
        + effective_potential(m=m, M=M, r=r, l=l)
        )

def effective_potential(
        m: float,
        M: float,
        r: float,
        l: float) -> float:
    return potential(m=m, M=M, r=r) + l**2 / (2. * m * r**2)

def theta_prime(m: float, l: float, r: float) -> float:
    return l / (m * r**2)

def r_prime(m: float, M: float, r: float, l: float, E: float) -> float:
    return np.sqrt(
        2. / m *(
            E - effective_potential(m=m, M=M, r=r, l=l)
        )
    )


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
        angular_momentum: float,
        m: float, M: float,
        total_energy: float) -> float:
    """https://de.wikipedia.org/wiki/Zweik%C3%B6rperproblem"""
    h = angular_momentum / m
    e = eccentricity(
        total_energy=total_energy,
        m=m, M=M, angular_momentum=angular_momentum)
    numerator = h**2 / (G*M)
    denominator = 1. + e * np.cos(theta - theta_0)
    return numerator / denominator
