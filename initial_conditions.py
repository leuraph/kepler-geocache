from scipy.constants import G


class InitialConditions:
    theta_0: float
    theta_prime_0: float
    r_0: float
    r_prime_0: float
    m: float
    M: float
    k: float
    v_0: float

    conserved_angular_momentum: float
    conserved_total_energy: float

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
        self.v_0 = theta_prime_0 * r_0

        self.conserved_angular_momentum = angular_momentum(
            m=m, r=r_0, theta_prime=theta_prime_0)
        self.conserved_total_energy = total_energy(
            m=m, M=M, r=r_0, r_prime=r_prime_0,
            l=self.conserved_angular_momentum)
    
    def get_conserved_angular_momentum(self) -> float:
        return self.conserved_angular_momentum
    
    def get_conserved_total_energy(self) -> float:
        return self.conserved_total_energy


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