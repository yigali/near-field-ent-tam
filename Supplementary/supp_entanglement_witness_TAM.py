import numpy as np
from numba import jit
import logging


logging.basicConfig(filename='dimensional_witness.log', level=logging.INFO, format='%(asctime)s - %(message)s')


@ jit(nopython=True)
def gamma(m, n, m_tag, n_tag, d):
    if (m - m_tag - n + n_tag) % d != 0:
        g = 0
    else:
        g = 1 / d
    return g


@ jit(nopython=True)
def fidelity1(rho):
    F1 = np.trace(rho) / rho.shape[0]
    return F1


@ jit(nopython=True)
def fidelity2(rhoR, rhoK):
    F2 = 0
    d = rhoR.shape[0]
    F2 += np.trace(rhoK)
    F2 -= 1 / d

    for m in range(d):
        for n in range(d):
            for m_tag in range(d):
                for n_tag in range(d):
                    if m == m_tag or m == n or n == n_tag or n_tag == m_tag:
                        continue
                    rhomn = rhoR[m, n] * rhoR[m_tag, n_tag]
                    sqrt_rho = np.sqrt(rhomn)
                    F2 -= gamma(m, n, m_tag, n_tag, d) * sqrt_rho
    return F2


def count_to_density_diagonal_matrix(count_matrix):
    # normalize the matrix which represent the main diagonal of the density matrix
    return count_matrix / count_matrix.sum()


def lambda_calc(rho, m):
    lambda_m = np.sqrt(rho[m, m] / np.trace(rho))
    return lambda_m


def c_lambda(rho):
    d = rho.shape[0]
    sum_lambda = 0
    c_lambda = 0
    for k in range(d):
        sum_lambda += lambda_calc(rho, k)
    for m in range(d):
        for n in range(d):
            c_lambda += lambda_calc(rho, m)*lambda_calc(rho, n)*rho[m,n]

    return ((d**2) / (sum_lambda**2)) * c_lambda


def DimensionalWitness(countK, countR):
    r = 0

    # count to density main diagonal data
    DdensityR = count_to_density_diagonal_matrix(countR)
    DdensityK = count_to_density_diagonal_matrix(countK)

    # c = c_lambda(DdensityR)
    # print(f'C_lambda = {c}')
    # DdensityK = DdensityK * c

    # Max entanglement dimension
    d_witness = DdensityR.shape[0]

    F1 = fidelity1(DdensityR)
    print(f'F1 is {F1}')
    F2 = fidelity2(DdensityR, DdensityK)
    print(f'F2 is {F2}')
    F = F1 + F2
    print(f'The estimate fidelity is {F}')
    while r / d_witness < F:
        r = r + 1
    print(f'The dimensional Witness is {r}')
    return r, F


def main():
    # test data
    # DdensityK = np.array([[1, 0],
    #                   [0, 1]])
    # DdensityR = np.array([[1, 0],
    #                   [0, 1]])

    # count matrices data
    cov_mat_K = np.fliplr(np.array([[16.0624, 177.5366],
                                 [181.3687, 14.9053]]))
    cov_mat_R = np.array([[156.1654, 9.56189],
                       [8.5684, 159.6432]])

    # Dimensional witness algorithm
    r, F = DimensionalWitness(cov_mat_K, cov_mat_R)
    logging.info(f'\n New run')
    logging.info(f'Dimensional Witness value: {r}')
    logging.info(f'Fidelity value: {F}')

if __name__ == "__main__":
    main()