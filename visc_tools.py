import numpy as np
from scipy.special import jv
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator as rgt

from math import cosh, sinh, cos, sin, radians, pi
import configparser

def get_constants():

    config = configparser.ConfigParser()
    config.read('params.cfg')
    viscoelastic_params = config["Viscoelastic Parameters"]

    vals = {'Nterms': int(viscoelastic_params["Nterms"]),
            'H1': float(viscoelastic_params["H1"]),
            'H2': float(viscoelastic_params["H2"]),
            'tR1': float(viscoelastic_params["tR1"]),
            'tR2': float(viscoelastic_params["tR2"]),
            'nu': float(viscoelastic_params["nu"]),
            'Gshear': float(viscoelastic_params["Gshear"]),
            'coseismic_slip': float(viscoelastic_params["coseismic_slip"]),
            't_eq_frac': float(viscoelastic_params["t_eq_frac"]),
            }

    vals['lam'] = 2 * vals['nu'] * vals['Gshear'] / (1 - (2 * vals['nu']))
    vals["bulk"] = vals['lam'] + 2 * vals["Gshear"] / 3
    # 
    vals['dH'] = vals['H1']-vals['H2']

    return vals


def momtensor_inverse(strike, dip, lam, mu):

    # Define unit vectors for strike and dip directions
    Vs = np.array([cos(strike), sin(strike), 0])
    Vd = np.array([-cos(dip) * sin(strike), cos(dip) * cos(strike), sin(dip)])

    # Fault-normal vector as cross product of dip and strike vectors
    Vnorm = np.cross(Vd, Vs)

    # Moment tensor for strike-slip component (rake = 0)
    M1 = np.zeros((3, 3))
    M1[0, 0] = (
        Vs[0] * Vnorm[0] * (lam + 2 * mu)
        + Vs[1] * Vnorm[1] * lam
        + Vs[2] * Vnorm[2] * lam
    )
    M1[1, 1] = (
        Vs[0] * Vnorm[0] * lam
        + Vs[1] * Vnorm[1] * (lam + 2 * mu)
        + Vs[2] * Vnorm[2] * lam
    )
    M1[2, 2] = (
        Vs[0] * Vnorm[0] * lam
        + Vs[1] * Vnorm[1] * lam
        + Vs[2] * Vnorm[2] * (lam + 2 * mu)
    )
    M1[0, 1] = M1[1, 0] = Vs[0] * Vnorm[1] * mu + Vs[1] * Vnorm[0] * mu
    M1[0, 2] = M1[2, 0] = Vs[0] * Vnorm[2] * mu + Vs[2] * Vnorm[0] * mu
    M1[1, 2] = M1[2, 1] = Vs[1] * Vnorm[2] * mu + Vs[2] * Vnorm[1] * mu

    # Moment tensor for dip-slip component (rake = 90)
    M2 = np.zeros((3, 3))
    M2[0, 0] = (
        Vd[0] * Vnorm[0] * (lam + 2 * mu)
        + Vd[1] * Vnorm[1] * lam
        + Vd[2] * Vnorm[2] * lam
    )
    M2[1, 1] = (
        Vd[0] * Vnorm[0] * lam
        + Vd[1] * Vnorm[1] * (lam + 2 * mu)
        + Vd[2] * Vnorm[2] * lam
    )
    M2[2, 2] = (
        Vd[0] * Vnorm[0] * lam
        + Vd[1] * Vnorm[1] * lam
        + Vd[2] * Vnorm[2] * (lam + 2 * mu)
    )
    M2[0, 1] = M2[1, 0] = Vd[0] * Vnorm[1] * mu + Vd[1] * Vnorm[0] * mu
    M2[0, 2] = M2[2, 0] = Vd[0] * Vnorm[2] * mu + Vd[2] * Vnorm[0] * mu
    M2[1, 2] = M2[2, 1] = Vd[1] * Vnorm[2] * mu + Vd[2] * Vnorm[1] * mu

    # Moment tensor for tensile component (fault opening)
    M3 = np.zeros((3, 3))
    M3[0, 0] = (
        Vnorm[0] * Vnorm[0] * (lam + 2 * mu)
        + Vnorm[1] * Vnorm[1] * lam
        + Vnorm[2] * Vnorm[2] * lam
    )
    M3[1, 1] = (
        Vnorm[0] * Vnorm[0] * lam
        + Vnorm[1] * Vnorm[1] * (lam + 2 * mu)
        + Vnorm[2] * Vnorm[2] * lam
    )
    M3[2, 2] = (
        Vnorm[0] * Vnorm[0] * lam
        + Vnorm[1] * Vnorm[1] * lam
        + Vnorm[2] * Vnorm[2] * (lam + 2 * mu)
    )
    M3[0, 1] = M3[1, 0] = Vnorm[0] * Vnorm[1] * mu + Vnorm[1] * Vnorm[0] * mu
    M3[0, 2] = M3[2, 0] = Vnorm[0] * Vnorm[2] * mu + Vnorm[2] * Vnorm[0] * mu
    M3[1, 2] = M3[2, 1] = Vnorm[1] * Vnorm[2] * mu + Vnorm[2] * Vnorm[1] * mu

    # Optional: Set small values close to zero
    M1[np.abs(M1) < 1e-13] = 0
    M2[np.abs(M2) < 1e-13] = 0
    M3[np.abs(M3) < 1e-13] = 0

    return M1, M2, M3


def get_prop(H, k, mu, lam, zs):
    g = lam + 2 * mu

    # Define the initial A matrix after Segall Eqn. 5.99
    # this is very similar to Eqn. 5.128 which has a few opposite signs?
    A4x4 = np.array(
        [
            [0, k, 1 / mu, 0],
            [-k * lam / g, 0, 0, 1 / g],
            [4 * k**2 * mu * (lam + mu) / g, 0, 0, k * lam / g],
            [0, 0, -k, 0],
        ]
    )

    # Compute constants for z=0 and z0=H
    z = 0
    z0 = H

    # simplify constants and propagator matrix after Segall eqn. 5.101
    C_pr = np.cosh(np.abs(k) * (z - z0))
    S_pr = np.sinh(np.abs(k) * (z - z0))

    # via the Cayley-Hamilton theorem, A can be expressed in terms of its own characteristic polynomial (degree = 4)
    # thus: e^(A*(z-z0)) = C_3 * A^3 + C_2 * A^2 + C_1 * A + C_0 * I
    # https://en.wikipedia.org/wiki/Cayley–Hamilton_theorem#Matrix_functions
    # has an example showing this for a 2x2 matrix
    C3 = -(S_pr - k * (z - z0) * C_pr) / (2 * k**3)
    C2 = k * (z - z0) * S_pr / (2 * k**2)
    C1 = (3 * S_pr - k * (z - z0) * C_pr) / (2 * k)
    C0 = (2 * C_pr - k * (z - z0) * S_pr) / 2

    # This is equivalent to Segall eqn. 5.100 where lambda = mu, but more generalzied
    # Propagator matrix, P4x4, is given as matrix exponential of A, that is P4x4 = e^A(z-z0)
    # while exp(A(z-z0)) is equal to the power series: n=0->n=oo sum: A*(z-z0)^n / n!
    P4x4 = (
        C3 * np.linalg.matrix_power(A4x4, 3)
        + C2 * np.linalg.matrix_power(A4x4, 2)
        + C1 * A4x4
        + C0 * np.eye(4)
    )

    # G matrix with small value for rg (regularization for numerical stability?)
    rg = 3 * 10**-3
    G = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, -rg, 0, 1]])

    # Segall eqn. 5.61

    P2x2 = np.array([[C_pr, (1 / (mu * abs(k))) * S_pr], [mu * abs(k) * S_pr, C_pr]])

    # Compute halfspace propagators
    halfspaceP4x4 = G @ P4x4
    halfspaceP2x2 = P2x2

    # Compute constants for z=0 and z0=zs for source propagator
    # simplify constants and propagator matrix after Segall eqn. 5.101
    C_pr = np.cosh(k * (z - zs))
    S_pr = np.sinh(k * (z - zs))

    C3 = -(S_pr - k * (z - zs) * C_pr) / (2 * k**3)
    C2 = k * (z - zs) * S_pr / (2 * k**2)
    C1 = (3 * S_pr - k * (z - zs) * C_pr) / (2 * k)
    C0 = (2 * C_pr - k * (z - zs) * S_pr) / 2

    # Compute P4x4zs and P2x2zs matrices
    P4x4zs = (
        C3 * np.linalg.matrix_power(A4x4, 3)
        + C2 * np.linalg.matrix_power(A4x4, 2)
        + C1 * A4x4
        + C0 * np.eye(4)
    )
    P2x2zs = np.array([[C_pr, (1 / (mu * abs(k))) * S_pr], [mu * abs(k) * S_pr, C_pr]])

    # Compute source propagators
    sourceP4x4 = G @ P4x4zs
    sourceP2x2 = P2x2zs

    return sourceP4x4, sourceP2x2, halfspaceP4x4, halfspaceP2x2


def get_denom_2visc(P0, d1, kj, B1, B2, K):

    P011, P012, P013, P014 = P0[0, :]
    P021, P022, P023, P024 = P0[1, :]
    P031, P032, P033, P034 = P0[2, :]
    P041, P042, P043, P044 = P0[3, :]

    # Calculating hyperbolic functions
    C1 = cosh(kj * d1)
    S1 = sinh(kj * d1)

    D1 = (
        64
        * B1**4
        * B2
        * K**3
        * (P032 * P041 - P031 * P042)
        * (
            -2 * B1 * B2 * C1 * S1
            + B2**2 * (C1**2 * (1 + d1**2 * kj**2) - d1**2 * kj**2 * S1**2)
            + B1**2 * (-(C1**2) * d1**2 * kj**2 + (1 + d1**2 * kj**2) * S1**2)
        )
    )

    D2 = (
        64
        * B1**3
        * K**2
        * (
            (-(B1**3))
            * (P032 * P041 - P031 * P042)
            * (C1**2 * d1**2 * kj**2 - (1 + d1**2 * kj**2) * S1**2)
            + B1
            * B2**2
            * (
                C1**2
                * (
                    (-(2 + K)) * P031 * P042
                    - d1**2 * kj**2 * (2 + K) * P031 * P042
                    + P032 * ((1 + d1**2 * kj**2) * (2 + K) * P041 - kj * K * P043)
                    + kj * K * ((-P034) * P041 + P033 * P042 + P031 * P044)
                )
                - C1 * (7 + 3 * K) * (P032 * P041 - P031 * P042) * S1
                + kj
                * (
                    (-(d1**2)) * kj * (2 + K) * (P032 * P041 - P031 * P042)
                    + K * ((-P034) * P041 + P033 * P042 - P032 * P043 + P031 * P044)
                )
                * S1**2
            )
            + B1**2
            * B2
            * (
                C1**2
                * d1
                * kj**2
                * (
                    K * (P034 * P041 + P033 * P042 - P032 * P043 - P031 * P044)
                    + d1
                    * (
                        (3 + 2 * K) * P031 * P042
                        + kj * K * (P033 * P041 - P034 * P042 - P031 * P043)
                        + P032 * (-(3 + 2 * K) * P041 + kj * K * P044)
                    )
                )
                + C1
                * (
                    (3 + K) * P031 * P042
                    + P032 * (-(3 + K) * P041 + kj * K * P043)
                    + kj * K * (P034 * P041 - P033 * P042 - P031 * P044)
                )
                * S1
                + (
                    -2 * (2 + K) * P031 * P042
                    - d1**2
                    * kj**2
                    * (
                        (3 + 2 * K) * P031 * P042
                        + kj * K * (P033 * P041 - P034 * P042 - P031 * P043)
                    )
                    - d1 * kj**2 * K * (P034 * P041 + P033 * P042 - P031 * P044)
                    + P032
                    * (
                        (2 * (2 + K) + d1**2 * kj**2 * (3 + 2 * K)) * P041
                        + d1 * kj**2 * K * (P043 - d1 * kj * P044)
                    )
                )
                * S1**2
            )
            + B2**3
            * (
                C1**2
                * (
                    (-(3 + K)) * P031 * P042
                    - d1**2
                    * kj**2
                    * (
                        (2 + K) * P031 * P042
                        + kj * K * (P033 * P041 - P034 * P042 - P031 * P043)
                    )
                    - d1 * kj**2 * K * (P034 * P041 + P033 * P042 - P031 * P044)
                    + P032
                    * (
                        (3 + K + d1**2 * kj**2 * (2 + K)) * P041
                        + d1 * kj**2 * K * (P043 - d1 * kj * P044)
                    )
                )
                + C1
                * kj
                * K
                * (P034 * P041 - P033 * P042 + P032 * P043 - P031 * P044)
                * S1
                + d1
                * kj**2
                * (
                    K * (P034 * P041 + P033 * P042 - P032 * P043 - P031 * P044)
                    + d1
                    * (
                        (2 + K) * P031 * P042
                        + kj * K * (P033 * P041 - P034 * P042 - P031 * P043)
                        + P032 * (-(2 + K) * P041 + kj * K * P044)
                    )
                )
                * S1**2
            )
        )
    )

    D3 = (
        -16
        * B1**2
        * K
        * (
            -2
            * B1**3
            * (
                C1**2
                * d1
                * kj**2
                * (
                    2 * K * (P034 * P041 + P033 * P042 - P032 * P043 - P031 * P044)
                    + d1
                    * (
                        (4 + 5 * K) * P031 * P042
                        + 2 * kj * K * (P033 * P041 - P034 * P042 - P031 * P043)
                        + P032 * (-(4 + 5 * K) * P041 + 2 * kj * K * P044)
                    )
                )
                + C1
                * K
                * (
                    -3 * P032 * P041
                    + 2 * kj * P034 * P041
                    + 3 * P031 * P042
                    - 2 * kj * P033 * P042
                    + 2 * kj * P032 * P043
                    - 2 * kj * P031 * P044
                )
                * S1
                + (
                    -4 * (2 + K) * P031 * P042
                    - d1**2
                    * kj**2
                    * (
                        (4 + 5 * K) * P031 * P042
                        + 2 * kj * K * (P033 * P041 - P034 * P042 - P031 * P043)
                    )
                    - 2 * d1 * kj**2 * K * (P034 * P041 + P033 * P042 - P031 * P044)
                    + P032
                    * (
                        (4 * (2 + K) + d1**2 * kj**2 * (4 + 5 * K)) * P041
                        + 2 * d1 * kj**2 * K * (P043 - d1 * kj * P044)
                    )
                )
                * S1**2
            )
            - 2
            * B1
            * B2**2
            * (
                2
                * C1**2
                * (
                    -(d1**2)
                    * kj**2
                    * (2 + K)
                    * (
                        (2 + K) * P031 * P042
                        + kj * K * (P033 * P041 - P034 * P042 - P031 * P043)
                    )
                    - d1
                    * kj**2
                    * K
                    * (2 + K)
                    * (P034 * P041 + P033 * P042 - P031 * P044)
                    + P032
                    * (
                        (2 + K) * (3 + K + d1**2 * kj**2 * (2 + K)) * P041
                        - kj
                        * K
                        * (
                            (3 + K - d1 * kj * (2 + K)) * P043
                            + d1**2 * kj**2 * (2 + K) * P044
                        )
                    )
                    - (3 + K)
                    * (
                        (2 + K) * P031 * P042
                        + kj * K * (P034 * P041 - P033 * P042 - P031 * P044)
                    )
                )
                + C1
                * (
                    3 * (4 + 7 * K + K**2) * P031 * P042
                    + 4 * kj**2 * K**2 * (-P034 * P043 + P033 * P044)
                    + P032
                    * (
                        -3 * (4 + 7 * K + K**2) * P041
                        + 2 * kj * K * ((2 + K) * P043 + P044)
                    )
                    - 2
                    * kj
                    * K
                    * (
                        P034 * (-(2 + K) * P041 + P042)
                        + P033 * (-P041 + (2 + K) * P042)
                        + P031 * (P043 + (2 + K) * P044)
                    )
                )
                * S1
                + 2
                * kj
                * (
                    d1
                    * kj
                    * K
                    * (2 + K)
                    * (P034 * P041 + P033 * P042 - P032 * P043 - P031 * P044)
                    - K
                    * (3 + K)
                    * (P034 * P041 - P033 * P042 + P032 * P043 - P031 * P044)
                    + d1**2
                    * kj
                    * (2 + K)
                    * (
                        (2 + K) * P031 * P042
                        + kj * K * (P033 * P041 - P034 * P042 - P031 * P043)
                        + P032 * (-(2 + K) * P041 + kj * K * P044)
                    )
                )
                * S1**2
            )
            + B2**3
            * (
                C1**2
                * (
                    (9 + 12 * K + K**2) * P031 * P042
                    + 2
                    * d1
                    * kj**2
                    * K
                    * (5 + K)
                    * (P034 * P041 + P033 * P042 - P031 * P044)
                    - P032
                    * (
                        (9 + 12 * K + K**2 + d1**2 * kj**2 * (4 + 8 * K + K**2)) * P041
                        + 2
                        * d1
                        * kj**2
                        * K
                        * ((5 + K) * P043 - d1 * kj * (4 + K) * P044)
                    )
                    + d1**2
                    * kj**2
                    * (
                        (4 + 8 * K + K**2) * P031 * P042
                        + 2
                        * kj
                        * K
                        * (4 + K)
                        * (P033 * P041 - P034 * P042 - P031 * P043)
                        + 4 * kj**2 * K**2 * (P034 * P043 - P033 * P044)
                    )
                )
                - 2
                * C1
                * kj
                * K
                * (5 + K)
                * (P034 * P041 - P033 * P042 + P032 * P043 - P031 * P044)
                * S1
                + (
                    -P031 * P042
                    + 2 * kj * K * (P033 * P041 - P034 * P042 - P031 * P043)
                    - 2
                    * d1**2
                    * kj**3
                    * K
                    * (4 + K)
                    * (P033 * P041 - P034 * P042 - P031 * P043)
                    + 4 * d1**2 * kj**4 * K**2 * (-P034 * P043 + P033 * P044)
                    + P032
                    * (
                        (1 + d1**2 * kj**2 * (4 + 8 * K + K**2)) * P041
                        + 2
                        * kj
                        * K
                        * (
                            d1 * kj * (5 + K) * P043
                            + P044
                            - d1**2 * kj**2 * (4 + K) * P044
                        )
                    )
                    - kj**2
                    * (
                        d1**2 * (4 + 8 * K + K**2) * P031 * P042
                        + 2
                        * d1
                        * K
                        * (5 + K)
                        * (P034 * P041 + P033 * P042 - P031 * P044)
                        + 4 * K**2 * (P034 * P043 - P033 * P044)
                    )
                )
                * S1**2
            )
            + (
                B1**2
                * B2
                * (
                    C1**2
                    * (
                        -2
                        * d1
                        * kj**2
                        * K
                        * (7 + 3 * K)
                        * (P034 * P041 + P033 * P042 - P032 * P043 - P031 * P044)
                        + d1**2
                        * kj**2
                        * (
                            -(12 + 14 * K + 5 * K**2) * P031 * P042
                            - 6
                            * kj
                            * K
                            * (2 + K)
                            * (P033 * P041 - P034 * P042 - P031 * P043)
                            + P032
                            * (
                                (12 + 14 * K + 5 * K**2) * P041
                                - 6 * kj * K * (2 + K) * P044
                            )
                            + 4 * kj**2 * K**2 * (-P034 * P043 + P033 * P044)
                        )
                        + K
                        * (
                            (8 + K) * P031 * P042
                            + kj**2 * K * (-4 * P034 * P043 + 4 * P033 * P044)
                            + P032
                            * (-(8 + K) * P041 + 2 * kj * ((3 + K) * P043 + P044))
                            + 2
                            * kj
                            * (
                                P034 * ((3 + K) * P041 - P042)
                                + P033 * (P041 - (3 + K) * P042)
                                - P031 * (P043 + (3 + K) * P044)
                            )
                        )
                    )
                    + 2
                    * C1
                    * (7 + 3 * K)
                    * (
                        -(3 + K) * P031 * P042
                        + P032 * ((3 + K) * P041 - kj * K * P043)
                        + kj * K * (-P034 * P041 + P033 * P042 + P031 * P044)
                    )
                    * S1
                    + (
                        2
                        * d1
                        * kj**2
                        * K
                        * (7 + 3 * K)
                        * (P034 * P041 + P033 * P042 - P031 * P044)
                        - P032
                        * (
                            (
                                18
                                + 32 * K
                                + 6 * K**2
                                + d1**2 * kj**2 * (12 + 14 * K + 5 * K**2)
                            )
                            * P041
                            - 2
                            * kj
                            * K
                            * (
                                (3 + K - d1 * kj * (7 + 3 * K)) * P043
                                + 3 * d1**2 * kj**2 * (2 + K) * P044
                            )
                        )
                        + 2
                        * (
                            (9 + 16 * K + 3 * K**2) * P031 * P042
                            + kj
                            * K
                            * (3 + K)
                            * (P034 * P041 - P033 * P042 - P031 * P044)
                        )
                        + d1**2
                        * kj**2
                        * (
                            (12 + 14 * K + 5 * K**2) * P031 * P042
                            + 6
                            * kj
                            * K
                            * (2 + K)
                            * (P033 * P041 - P034 * P042 - P031 * P043)
                            + 4 * kj**2 * K**2 * (P034 * P043 - P033 * P044)
                        )
                    )
                    * S1**2
                )
            )
        )
    )

    D4 = (
        -16
        * B1
        * (
            B1**2
            * B2
            * (
                (
                    C1**2
                    * (
                        (-d1)
                        * kj**2
                        * K
                        * (16 + 11 * K + 2 * K**2)
                        * (P034 * P041 + P033 * P042 - P032 * P043 - P031 * P044)
                        + d1**2
                        * kj**2
                        * (
                            (-(4 + 4 * K + 5 * K**2 + K**3)) * P031 * P042
                            - 2
                            * kj
                            * K
                            * (6 + 4 * K + K**2)
                            * (P033 * P041 - P034 * P042 - P031 * P043)
                            + P032
                            * (
                                (4 + 4 * K + 5 * K**2 + K**3) * P041
                                - 2 * kj * K * (6 + 4 * K + K**2) * P044
                            )
                            - 4 * kj**2 * K**2 * (3 + K) * (P034 * P043 - P033 * P044)
                        )
                        - K
                        * (3 + K)
                        * (
                            (-(8 + K)) * P031 * P042
                            + 4 * kj**2 * K * (P034 * P043 - P033 * P044)
                            + P032 * ((8 + K) * P041 - 2 * kj * ((3 + K) * P043 + P044))
                            + 2
                            * kj
                            * (
                                P034 * (-(3 + K) * P041 + P042)
                                + P033 * (-P041 + (3 + K) * P042)
                                + P031 * (P043 + (3 + K) * P044)
                            )
                        )
                    )
                )
                + (
                    C1
                    * (
                        -3 * (12 + 25 * K + 10 * K**2 + K**3) * P031 * P042
                        + 4 * kj**2 * K**2 * (3 + K) * (P034 * P043 - P033 * P044)
                        + P032
                        * (
                            3 * (12 + 25 * K + 10 * K**2 + K**3) * P041
                            - kj
                            * K
                            * ((14 + 29 * K + 4 * K**2) * P043 + 2 * (3 + K) * P044)
                        )
                        + kj
                        * K
                        * (
                            P034
                            * (-(14 + 29 * K + 4 * K**2) * P041 + 2 * (3 + K) * P042)
                            + P033
                            * (-2 * (3 + K) * P041 + (14 + 29 * K + 4 * K**2) * P042)
                            + P031
                            * (2 * (3 + K) * P043 + (14 + 29 * K + 4 * K**2) * P044)
                        )
                    )
                    * S1
                )
                + (
                    (4 + 35 * K + 24 * K**2 + 2 * K**3) * P031 * P042
                    + 2
                    * d1**2
                    * kj**3
                    * K
                    * (6 + 4 * K + K**2)
                    * (P033 * P041 - P034 * P042 - P031 * P043)
                    + 4 * d1**2 * kj**4 * K**2 * (3 + K) * (P034 * P043 - P033 * P044)
                    - kj
                    * K
                    * (
                        (-P034) * (2 * (3 + K) ** 2 * P041 - (-2 + K) * P042)
                        + P033 * ((-(-2 + K)) * P041 + 2 * (3 + K) ** 2 * P042)
                        + P031 * ((-2 + K) * P043 + 2 * (3 + K) ** 2 * P044)
                    )
                    - P032
                    * (
                        (
                            4
                            + 35 * K
                            + 24 * K**2
                            + 2 * K**3
                            + d1**2 * kj**2 * (4 + 4 * K + 5 * K**2 + K**3)
                        )
                        * P041
                        - kj
                        * K
                        * (
                            (2 * (3 + K) ** 2 - d1 * kj * (16 + 11 * K + 2 * K**2))
                            * P043
                            + (-2 + K + 2 * d1**2 * kj**2 * (6 + 4 * K + K**2)) * P044
                        )
                    )
                    + d1
                    * kj**2
                    * (
                        4 * d1 * P031 * P042
                        + K**2
                        * (
                            11 * P034 * P041
                            + 5 * d1 * P031 * P042
                            + 11 * P033 * P042
                            - 11 * P031 * P044
                        )
                        + 4
                        * K
                        * (
                            4 * P034 * P041
                            + d1 * P031 * P042
                            + 4 * P033 * P042
                            - 4 * P031 * P044
                        )
                        + K**3
                        * (
                            2 * P034 * P041
                            + d1 * P031 * P042
                            + 2 * P033 * P042
                            - 2 * P031 * P044
                        )
                    )
                )
                * S1**2
            )
            + (
                B1**3
                * (
                    C1**2
                    * (
                        -2
                        * d1
                        * kj**2
                        * K
                        * (5 + 4 * K)
                        * (P034 * P041 + P033 * P042 - P032 * P043 - P031 * P044)
                        + d1**2
                        * kj**2
                        * (
                            -(4 + 20 * K + 7 * K**2) * P031 * P042
                            - 8
                            * kj
                            * K
                            * (1 + K)
                            * (P033 * P041 - P034 * P042 - P031 * P043)
                            + P032
                            * (
                                (4 + 20 * K + 7 * K**2) * P041
                                - 8 * kj * K * (1 + K) * P044
                            )
                            + 4 * kj**2 * K**2 * (-P034 * P043 + P033 * P044)
                        )
                        + K**2
                        * (
                            2 * P031 * P042
                            + kj
                            * (
                                P033 * P041
                                + 3 * P034 * P041
                                - 3 * P033 * P042
                                - P034 * P042
                                - P031 * P043
                                - 3 * P031 * P044
                            )
                            + kj**2 * (-4 * P034 * P043 + 4 * P033 * P044)
                            + P032 * (-2 * P041 + kj * (3 * P043 + P044))
                        )
                    )
                    + C1
                    * K
                    * (7 + 3 * K)
                    * (
                        -3 * P031 * P042
                        + P032 * (3 * P041 - 2 * kj * P043)
                        + 2 * kj * (-P034 * P041 + P033 * P042 + P031 * P044)
                    )
                    * S1
                    + (
                        (16 + 33 * K + 6 * K**2) * P031 * P042
                        + 8
                        * d1**2
                        * kj**3
                        * K
                        * (1 + K)
                        * (P033 * P041 - P034 * P042 - P031 * P043)
                        + 4 * d1**2 * kj**4 * K**2 * (P034 * P043 - P033 * P044)
                        - P032
                        * (
                            (
                                16
                                + 33 * K
                                + 6 * K**2
                                + d1**2 * kj**2 * (4 + 20 * K + 7 * K**2)
                            )
                            * P041
                            + kj
                            * K
                            * (
                                -3 * K * P043
                                + 2 * d1 * kj * (5 + 4 * K) * P043
                                - 2 * P044
                                + K * P044
                                - 8 * d1**2 * kj**2 * (1 + K) * P044
                            )
                        )
                        + kj
                        * K
                        * (
                            -P033 * ((-2 + K) * P041 + 3 * K * P042)
                            - 2 * (P034 * P042 + P031 * P043)
                            + K
                            * (
                                3 * P034 * P041
                                + P034 * P042
                                + P031 * P043
                                - 3 * P031 * P044
                            )
                        )
                        + d1
                        * kj**2
                        * (
                            4 * d1 * P031 * P042
                            + K**2
                            * (
                                8 * P034 * P041
                                + 7 * d1 * P031 * P042
                                + 8 * P033 * P042
                                - 8 * P031 * P044
                            )
                            + 10
                            * K
                            * (
                                P034 * P041
                                + 2 * d1 * P031 * P042
                                + P033 * P042
                                - P031 * P044
                            )
                        )
                    )
                    * S1**2
                )
            )
            + (
                B2**3
                * K
                * (
                    C1**2
                    * (
                        3 * (3 + K) * P031 * P042
                        + d1
                        * kj**2
                        * (6 + 5 * K)
                        * (P034 * P041 + P033 * P042 - P031 * P044)
                        - P032
                        * (
                            (2 * d1**2 * kj**2 * (2 + K) + 3 * (3 + K)) * P041
                            + d1
                            * kj**2
                            * ((6 + 5 * K) * P043 - 4 * d1 * kj * (1 + K) * P044)
                        )
                        + 2
                        * d1**2
                        * kj**2
                        * (
                            (2 + K) * P031 * P042
                            + 2
                            * kj
                            * (1 + K)
                            * (P033 * P041 - P034 * P042 - P031 * P043)
                            + 4 * kj**2 * K * (P034 * P043 - P033 * P044)
                        )
                    )
                    - C1
                    * kj
                    * (6 + 5 * K)
                    * (P034 * P041 - P033 * P042 + P032 * P043 - P031 * P044)
                    * S1
                    + (
                        -P031 * P042
                        - 4
                        * d1**2
                        * kj**3
                        * (1 + K)
                        * (P033 * P041 - P034 * P042 - P031 * P043)
                        + kj * (2 + K) * (P033 * P041 - P034 * P042 - P031 * P043)
                        + 8 * d1**2 * kj**4 * K * (-P034 * P043 + P033 * P044)
                        + P032
                        * (
                            (1 + 2 * d1**2 * kj**2 * (2 + K)) * P041
                            + kj
                            * (
                                d1 * kj * (6 + 5 * K) * P043
                                - 4 * d1**2 * kj**2 * (1 + K) * P044
                                + (2 + K) * P044
                            )
                        )
                        - kj**2
                        * (
                            2 * d1**2 * (2 + K) * P031 * P042
                            + d1
                            * (6 + 5 * K)
                            * (P034 * P041 + P033 * P042 - P031 * P044)
                            + 8 * K * (P034 * P043 - P033 * P044)
                        )
                    )
                    * S1**2
                )
            )
            + B1
            * B2**2
            * (
                (
                    C1**2
                    * (
                        2
                        * d1
                        * kj**2
                        * K
                        * (10 + 7 * K + K**2)
                        * (P034 * P041 + P033 * P042 - P031 * P044)
                        - P032
                        * (
                            (2 + K)
                            * (9 + 12 * K + K**2 + d1**2 * kj**2 * (4 + 8 * K + K**2))
                            * P041
                            - kj
                            * K
                            * (
                                (9 + 12 * K + K**2 - 2 * d1 * kj * (10 + 7 * K + K**2))
                                * P043
                                + 2 * d1**2 * kj**2 * (8 + 6 * K + K**2) * P044
                            )
                        )
                        + (9 + 12 * K + K**2)
                        * (
                            (2 + K) * P031 * P042
                            + kj * K * (P034 * P041 - P033 * P042 - P031 * P044)
                        )
                        - d1**2
                        * kj**2
                        * (2 + K)
                        * (
                            -(4 + 8 * K + K**2) * P031 * P042
                            - 2
                            * kj
                            * K
                            * (4 + K)
                            * (P033 * P041 - P034 * P042 - P031 * P043)
                            + 4 * kj**2 * K**2 * ((-P034) * P043 + P033 * P044)
                        )
                    )
                )
                + (
                    C1
                    * K
                    * (
                        (-(36 + 21 * K + K**2)) * P031 * P042
                        + 4 * kj**2 * K * (5 + K) * (P034 * P043 - P033 * P044)
                        + P032
                        * (
                            (36 + 21 * K + K**2) * P041
                            - 2 * kj * ((10 + 7 * K + K**2) * P043 + (3 + 2 * K) * P044)
                        )
                        + 2
                        * kj
                        * (
                            -P034 * ((10 + 7 * K + K**2) * P041 - (3 + 2 * K) * P042)
                            + P033 * (-(3 + 2 * K) * P041 + (10 + 7 * K + K**2) * P042)
                            + P031 * ((3 + 2 * K) * P043 + (10 + 7 * K + K**2) * P044)
                        )
                    )
                    * S1
                )
                + (
                    (-(2 + K)) * P031 * P042
                    - 2
                    * d1**2
                    * kj**3
                    * K
                    * (8 + 6 * K + K**2)
                    * (P033 * P041 - P034 * P042 - P031 * P043)
                    - 4 * d1**2 * kj**4 * K**2 * (2 + K) * (P034 * P043 - P033 * P044)
                    - kj
                    * K
                    * (
                        -P034 * ((9 + 12 * K + K**2) * P041 - 2 * (2 + K) * P042)
                        + P033 * (-2 * (2 + K) * P041 + (9 + 12 * K + K**2) * P042)
                        + P031 * (2 * (2 + K) * P043 + (9 + 12 * K + K**2) * P044)
                    )
                    + P032
                    * (
                        (2 + K) * (1 + d1**2 * kj**2 * (4 + 8 * K + K**2)) * P041
                        + kj
                        * K
                        * (
                            (9 + 12 * K + K**2 + 2 * d1 * kj * (10 + 7 * K + K**2))
                            * P043
                            - 2 * (2 + K) * (-1 + d1**2 * kj**2 * (4 + K)) * P044
                        )
                    )
                    - kj**2
                    * (2 + K)
                    * (
                        d1**2 * (4 + 8 * K + K**2) * P031 * P042
                        + 2
                        * d1
                        * K
                        * (5 + K)
                        * (P034 * P041 + P033 * P042 - P031 * P044)
                        + 4 * K**2 * (P034 * P043 - P033 * P044)
                    )
                )
                * S1**2
            )
        )
    )

    D5 = -4 * (
        (
            B2**3
            * K
            * (
                C1**2
                * (
                    9 * P031 * P042
                    + 12 * d1 * kj**2 * (P034 * P041 + P033 * P042 - P031 * P044)
                    - P032
                    * (
                        (9 + 4 * d1**2 * kj**2) * P041
                        + 4 * d1 * kj**2 * (3 * P043 - 2 * d1 * kj * P044)
                    )
                    - 4
                    * d1**2
                    * kj**2
                    * (
                        (-P031) * P042
                        + 2 * kj * ((-P033) * P041 + P034 * P042 + P031 * P043)
                        + kj**2 * (-4 * P034 * P043 + 4 * P033 * P044)
                    )
                )
                - 12
                * C1
                * kj
                * (P034 * P041 - P033 * P042 + P032 * P043 - P031 * P044)
                * S1
                + (
                    (-P031) * P042
                    + 4 * kj * (P033 * P041 - P034 * P042 - P031 * P043)
                    + 8 * d1**2 * kj**3 * ((-P033) * P041 + P034 * P042 + P031 * P043)
                    - 16 * d1**2 * kj**4 * (P034 * P043 - P033 * P044)
                    + P032
                    * (
                        P041
                        + 4 * d1**2 * kj**2 * P041
                        + 4
                        * kj
                        * (3 * d1 * kj * P043 + P044 - 2 * d1**2 * kj**2 * P044)
                    )
                    - 4
                    * kj**2
                    * (
                        d1**2 * P031 * P042
                        + 4 * P034 * P043
                        - 4 * P033 * P044
                        + 3 * d1 * (P034 * P041 + P033 * P042 - P031 * P044)
                    )
                )
                * S1**2
            )
        )
        + (
            2
            * B1**3
            * (
                C1**2
                * (
                    -2
                    * d1
                    * kj**2
                    * (6 + 20 * K + 3 * K**2)
                    * (P034 * P041 + P033 * P042 - P032 * P043 - P031 * P044)
                    + d1**2
                    * kj**2
                    * (
                        -(20 + 28 * K + 3 * K**2) * P031 * P042
                        - 2
                        * kj
                        * (4 + 16 * K + 3 * K**2)
                        * (P033 * P041 - P034 * P042 - P031 * P043)
                        + P032
                        * (
                            (20 + 28 * K + 3 * K**2) * P041
                            - 2 * kj * (4 + 16 * K + 3 * K**2) * P044
                        )
                        - 4 * kj**2 * K * (4 + 3 * K) * (P034 * P043 - P033 * P044)
                    )
                    + 2
                    * K
                    * (3 + K)
                    * (
                        2 * P031 * P042
                        + kj
                        * (
                            P033 * P041
                            + 3 * P034 * P041
                            - 3 * P033 * P042
                            - P034 * P042
                            - P031 * P043
                            - 3 * P031 * P044
                        )
                        + kj**2 * (-4 * P034 * P043 + 4 * P033 * P044)
                        + P032 * (-2 * P041 + kj * (3 * P043 + P044))
                    )
                )
                + C1
                * (
                    -9 * (4 + 7 * K + K**2) * P031 * P042
                    + 12 * kj**2 * K**2 * (P034 * P043 - P033 * P044)
                    + P032
                    * (
                        9 * (4 + 7 * K + K**2) * P041
                        - 2 * kj * ((12 + 22 * K + 5 * K**2) * P043 + 3 * K * P044)
                    )
                    + 2
                    * kj
                    * (
                        12 * ((-P034) * P041 + P033 * P042 + P031 * P044)
                        + 5 * K**2 * ((-P034) * P041 + P033 * P042 + P031 * P044)
                        + K
                        * (
                            -3 * P033 * P041
                            - 22 * P034 * P041
                            + 22 * P033 * P042
                            + 3 * P034 * P042
                            + 3 * P031 * P043
                            + 22 * P031 * P044
                        )
                    )
                )
                * S1
                + (
                    (68 + 51 * K + 4 * K**2) * P031 * P042
                    + 2
                    * d1**2
                    * kj**3
                    * (4 + 16 * K + 3 * K**2)
                    * (P033 * P041 - P034 * P042 - P031 * P043)
                    + 4 * d1**2 * kj**4 * K * (4 + 3 * K) * (P034 * P043 - P033 * P044)
                    - P032
                    * (
                        (
                            68
                            + 51 * K
                            + 4 * K**2
                            + d1**2 * kj**2 * (20 + 28 * K + 3 * K**2)
                        )
                        * P041
                        - 2
                        * kj
                        * (
                            9 * K * P043
                            - d1 * kj * (6 + 20 * K + 3 * K**2) * P043
                            + K**2 * (3 * P043 - P044)
                            + 4 * P044
                            + d1**2 * kj**2 * (4 + 16 * K + 3 * K**2) * P044
                        )
                    )
                    + d1
                    * kj**2
                    * (
                        2 * (6 + 20 * K + 3 * K**2) * P034 * P041
                        + d1 * (20 + 28 * K + 3 * K**2) * P031 * P042
                        + 2 * (6 + 20 * K + 3 * K**2) * (P033 * P042 - P031 * P044)
                    )
                    - 2
                    * kj
                    * (
                        P033 * ((-4 + K**2) * P041 + 3 * K * (3 + K) * P042)
                        + 4 * (P034 * P042 + P031 * P043)
                        - K**2
                        * (
                            3 * P034 * P041
                            + P034 * P042
                            + P031 * P043
                            - 3 * P031 * P044
                        )
                        + K * (-9 * P034 * P041 + 9 * P031 * P044)
                    )
                )
                * S1**2
            )
        )
        + B1**2
        * B2
        * (
            (
                C1**2
                * (
                    (72 + 105 * K + 20 * K**2 + K**3) * P031 * P042
                    - 4
                    * d1**2
                    * kj**3
                    * (4 - 2 * K + K**2)
                    * (P033 * P041 - P034 * P042 - P031 * P043)
                    - 8 * d1**2 * kj**4 * K * (6 + K) * (P034 * P043 - P033 * P044)
                    - 2
                    * kj
                    * (9 + 12 * K + K**2)
                    * (
                        P034 * (-(3 + K) * P041 + P042)
                        + P033 * (-P041 + (3 + K) * P042)
                        + P031 * (P043 + (3 + K) * P044)
                    )
                    + P032
                    * (
                        (
                            -72
                            - 105 * K
                            - 20 * K**2
                            - K**3
                            + 2 * d1**2 * kj**2 * (-4 - 10 * K + K**2)
                        )
                        * P041
                        + 2
                        * kj
                        * (
                            (
                                27
                                + 45 * K
                                + 15 * K**2
                                + K**3
                                + 2 * d1 * kj * (6 - K + 2 * K**2)
                            )
                            * P043
                            + (
                                9
                                + 12 * K
                                + K**2
                                - 2 * d1**2 * kj**2 * (4 - 2 * K + K**2)
                            )
                            * P044
                        )
                    )
                    + kj**2
                    * (
                        -2 * d1**2 * (-4 - 10 * K + K**2) * P031 * P042
                        - 4
                        * d1
                        * (6 - K + 2 * K**2)
                        * (P034 * P041 + P033 * P042 - P031 * P044)
                        - 4 * K * (9 + 12 * K + K**2) * (P034 * P043 - P033 * P044)
                    )
                )
            )
            + (
                2
                * C1
                * (
                    (-(108 + 99 * K + 24 * K**2 + K**3)) * P031 * P042
                    + 4 * kj**2 * K * (15 + 8 * K + K**2) * (P034 * P043 - P033 * P044)
                    + P032
                    * (
                        (108 + 99 * K + 24 * K**2 + K**3) * P041
                        - 2
                        * kj
                        * (
                            (3 + 40 * K + 17 * K**2 + K**3) * P043
                            + (9 + 9 * K + 2 * K**2) * P044
                        )
                    )
                    + 2
                    * kj
                    * (
                        (-P034)
                        * (
                            (3 + 40 * K + 17 * K**2 + K**3) * P041
                            - (9 + 9 * K + 2 * K**2) * P042
                        )
                        + P033
                        * (
                            -(9 + 9 * K + 2 * K**2) * P041
                            + (3 + 40 * K + 17 * K**2 + K**3) * P042
                        )
                        + P031
                        * (
                            (9 + 9 * K + 2 * K**2) * P043
                            + (3 + 40 * K + 17 * K**2 + K**3) * P044
                        )
                    )
                )
            )
            * S1
            + (
                (16 + 101 * K + 32 * K**2 + K**3) * P031 * P042
                + 4
                * d1**2
                * kj**3
                * (4 - 2 * K + K**2)
                * (P033 * P041 - P034 * P042 - P031 * P043)
                + 8 * d1**2 * kj**4 * K * (6 + K) * (P034 * P043 - P033 * P044)
                + 2
                * kj
                * (
                    P034
                    * (
                        (27 + 45 * K + 15 * K**2 + K**3) * P041
                        + (7 - 8 * K - 3 * K**2) * P042
                    )
                    + P033
                    * (
                        (-7 + 8 * K + 3 * K**2) * P041
                        - (27 + 45 * K + 15 * K**2 + K**3) * P042
                    )
                    - P031
                    * (
                        (-7 + 8 * K + 3 * K**2) * P043
                        + (27 + 45 * K + 15 * K**2 + K**3) * P044
                    )
                )
                + P032
                * (
                    (
                        -(
                            16
                            + 101 * K
                            + 32 * K**2
                            + K**3
                            + 2 * d1**2 * kj**2 * (-4 - 10 * K + K**2)
                        )
                    )
                    * P041
                    + 2
                    * kj
                    * (
                        (
                            27
                            + 45 * K
                            + 15 * K**2
                            + K**3
                            - 2 * d1 * kj * (6 - K + 2 * K**2)
                        )
                        * P043
                        + (
                            -7
                            + 8 * K
                            + 3 * K**2
                            + 2 * d1**2 * kj**2 * (4 - 2 * K + K**2)
                        )
                        * P044
                    )
                )
                + 2
                * kj**2
                * (
                    d1**2 * (-4 - 10 * K + K**2) * P031 * P042
                    + 2
                    * d1
                    * (6 - K + 2 * K**2)
                    * (P034 * P041 + P033 * P042 - P031 * P044)
                    - 2 * K * (1 + 8 * K + K**2) * (P034 * P043 - P033 * P044)
                )
            )
            * S1**2
        )
        - 2
        * B1
        * B2**2
        * (
            2
            * C1**2
            * (
                -d1
                * kj**2
                * (12 + 16 * K + 5 * K**2)
                * (P034 * P041 + P033 * P042 - P031 * P044)
                + P032
                * (
                    (2 + K) * (2 * d1**2 * kj**2 * (2 + K) + 3 * (3 + K)) * P041
                    - kj
                    * (
                        3 * K * (3 + K) * P043
                        - d1 * kj * (12 + 16 * K + 5 * K**2) * P043
                        + 4 * d1**2 * kj**2 * (2 + 3 * K + K**2) * P044
                    )
                )
                - 3
                * (3 + K)
                * (
                    (2 + K) * P031 * P042
                    + kj * K * (P034 * P041 - P033 * P042 - P031 * P044)
                )
                + 2
                * d1**2
                * kj**2
                * (2 + K)
                * (
                    -(2 + K) * P031 * P042
                    - 2 * kj * (1 + K) * (P033 * P041 - P034 * P042 - P031 * P043)
                    + kj**2 * K * (-4 * P034 * P043 + 4 * P033 * P044)
                )
            )
            + C1
            * (
                24 * kj * (P034 * P041 - P033 * P042 + P032 * P043 - P031 * P044)
                + 4
                * K
                * (
                    9 * P031 * P042
                    + P032 * (-9 * P041 + 8 * kj * P043 + 3 * kj * P044)
                    + kj
                    * (
                        3 * P033 * P041
                        + 8 * P034 * P041
                        - 8 * P033 * P042
                        - 3 * P034 * P042
                        - 3 * P031 * P043
                        - 8 * P031 * P044
                    )
                    + kj**2 * (-6 * P034 * P043 + 6 * P033 * P044)
                )
                + K**2
                * (
                    7 * P031 * P042
                    + 2
                    * kj
                    * (
                        P033 * P041
                        + 5 * P034 * P041
                        - 5 * P033 * P042
                        - P034 * P042
                        - P031 * P043
                        - 5 * P031 * P044
                    )
                    - 20 * kj**2 * (P034 * P043 - P033 * P044)
                    + P032 * (-7 * P041 + 2 * kj * (5 * P043 + P044))
                )
            )
            * S1
            - 2
            * (
                -(2 + K) * P031 * P042
                - 4
                * d1**2
                * kj**3
                * (2 + 3 * K + K**2)
                * (P033 * P041 - P034 * P042 - P031 * P043)
                - 8 * d1**2 * kj**4 * K * (2 + K) * (P034 * P043 - P033 * P044)
                + kj
                * (
                    P033 * ((2 + K) ** 2 * P041 - 3 * K * (3 + K) * P042)
                    - 4 * (P034 * P042 + P031 * P043)
                    + K
                    * (
                        9 * P034 * P041
                        - 4 * P034 * P042
                        - 4 * P031 * P043
                        - 9 * P031 * P044
                    )
                    + K**2
                    * (3 * P034 * P041 - P034 * P042 - P031 * P043 - 3 * P031 * P044)
                )
                - kj**2
                * (2 + K)
                * (
                    2 * d1**2 * (2 + K) * P031 * P042
                    + d1 * (6 + 5 * K) * (P034 * P041 + P033 * P042 - P031 * P044)
                    + 8 * K * (P034 * P043 - P033 * P044)
                )
                + P032
                * (
                    (2 + K) * (1 + 2 * d1**2 * kj**2 * (2 + K)) * P041
                    + kj
                    * (
                        d1 * kj * (12 + 16 * K + 5 * K**2) * P043
                        + 4 * P044
                        - 4 * d1**2 * kj**2 * (2 + 3 * K + K**2) * P044
                        + K**2 * (3 * P043 + P044)
                        + K * (9 * P043 + 4 * P044)
                    )
                )
            )
            * S1**2
        )
    )

    D6 = -4 * (
        B2**2
        * (
            C1**2
            * (
                12 * d1 * kj**2 * (2 + K) * (P034 * P041 + P033 * P042 - P031 * P044)
                - P032
                * (
                    (9 + 4 * d1**2 * kj**2) * (2 + K) * P041
                    + kj
                    * (
                        -9 * K * P043
                        + 12 * d1 * kj * (2 + K) * P043
                        - 8 * d1**2 * kj**2 * (2 + K) * P044
                    )
                )
                + 9
                * (
                    (2 + K) * P031 * P042
                    + kj * K * (P034 * P041 - P033 * P042 - P031 * P044)
                )
                - 4
                * d1**2
                * kj**2
                * (2 + K)
                * (
                    (-P031) * P042
                    + 2 * kj * ((-P033) * P041 + P034 * P042 + P031 * P043)
                    + kj**2 * (-4 * P034 * P043 + 4 * P033 * P044)
                )
            )
            + 6
            * C1
            * (
                4 * kj * ((-P034) * P041 + P033 * P042 - P032 * P043 + P031 * P044)
                + K
                * (
                    -2 * P031 * P042
                    + kj
                    * (
                        (-P033) * P041
                        - 2 * P034 * P041
                        + 2 * P033 * P042
                        + P034 * P042
                        + P031 * P043
                        + 2 * P031 * P044
                    )
                    + 4 * kj**2 * (P034 * P043 - P033 * P044)
                    + P032 * (2 * P041 - kj * (2 * P043 + P044))
                )
            )
            * S1
            + (
                (-(2 + K)) * P031 * P042
                - 8
                * d1**2
                * kj**3
                * (2 + K)
                * (P033 * P041 - P034 * P042 - P031 * P043)
                - 16 * d1**2 * kj**4 * (2 + K) * (P034 * P043 - P033 * P044)
                + P032
                * (
                    (1 + 4 * d1**2 * kj**2) * (2 + K) * P041
                    + kj
                    * (
                        9 * K * P043
                        + 12 * d1 * kj * (2 + K) * P043
                        + 8 * P044
                        + 4 * K * P044
                        - 8 * d1**2 * kj**2 * (2 + K) * P044
                    )
                )
                + kj
                * (
                    P033 * (4 * (2 + K) * P041 - 9 * K * P042)
                    - 8 * (P034 * P042 + P031 * P043)
                    + K
                    * (
                        9 * P034 * P041
                        - 4 * P034 * P042
                        - 4 * P031 * P043
                        - 9 * P031 * P044
                    )
                )
                - 4
                * kj**2
                * (2 + K)
                * (
                    d1**2 * P031 * P042
                    + 4 * P034 * P043
                    - 4 * P033 * P044
                    + 3 * d1 * (P034 * P041 + P033 * P042 - P031 * P044)
                )
            )
            * S1**2
        )
        + B1
        * B2
        * (
            C1**2
            * (
                3 * (24 + 11 * K + K**2) * P031 * P042
                + 16
                * d1**2
                * kj**3
                * (1 + K)
                * (P033 * P041 - P034 * P042 - P031 * P043)
                + 16 * d1**2 * kj**4 * (-1 + 2 * K) * (P034 * P043 - P033 * P044)
                - 6
                * kj
                * (3 + K)
                * (
                    P034 * (-(3 + K) * P041 + P042)
                    + P033 * (-P041 + (3 + K) * P042)
                    + P031 * (P043 + (3 + K) * P044)
                )
                + P032
                * (
                    -(4 * d1**2 * kj**2 * (5 + 2 * K) + 3 * (24 + 11 * K + K**2)) * P041
                    + 2
                    * kj
                    * (
                        -3 * (-((3 + K) ** 2) + d1 * kj * (4 + 3 * K)) * P043
                        + (8 * d1**2 * kj**2 * (1 + K) + 3 * (3 + K)) * P044
                    )
                )
                + 2
                * kj**2
                * (
                    2 * d1**2 * (5 + 2 * K) * P031 * P042
                    + 3 * d1 * (4 + 3 * K) * (P034 * P041 + P033 * P042 - P031 * P044)
                    - 6 * K * (3 + K) * (P034 * P043 - P033 * P044)
                )
            )
            + C1
            * (
                (-(108 + 57 * K + 7 * K**2)) * P031 * P042
                + 4 * kj**2 * (18 + 21 * K + 5 * K**2) * (P034 * P043 - P033 * P044)
                + 2
                * kj
                * (
                    P033 * (-(18 + 9 * K + K**2) * P041 + 6 * (5 + 7 * K + K**2) * P042)
                    + P034
                    * (-6 * (5 + 7 * K + K**2) * P041 + (18 + 9 * K + K**2) * P042)
                    + P031
                    * ((18 + 9 * K + K**2) * P043 + 6 * (5 + 7 * K + K**2) * P044)
                )
                + P032
                * (
                    (108 + 57 * K + 7 * K**2) * P041
                    - 2
                    * kj
                    * (6 * (5 + 7 * K + K**2) * P043 + (18 + 9 * K + K**2) * P044)
                )
            )
            * S1
            + 2
            * (
                2 * (1 + 8 * K + K**2) * P031 * P042
                - 8
                * d1**2
                * kj**3
                * (1 + K)
                * (P033 * P041 - P034 * P042 - P031 * P043)
                - 8 * d1**2 * kj**4 * (-1 + 2 * K) * (P034 * P043 - P033 * P044)
                + kj
                * (
                    P033 * ((1 + 8 * K + K**2) * P041 - 3 * (3 + K) ** 2 * P042)
                    + P034 * (3 * (3 + K) ** 2 * P041 - (1 + 8 * K + K**2) * P042)
                    - P031 * ((1 + 8 * K + K**2) * P043 + 3 * (3 + K) ** 2 * P044)
                )
                + P032
                * (
                    2 * (-1 - 8 * K - K**2 + d1**2 * kj**2 * (5 + 2 * K)) * P041
                    + kj
                    * (
                        3 * ((3 + K) ** 2 + d1 * kj * (4 + 3 * K)) * P043
                        + (1 + 8 * K + K**2 - 8 * d1**2 * kj**2 * (1 + K)) * P044
                    )
                )
                - kj**2
                * (
                    2 * d1**2 * (5 + 2 * K) * P031 * P042
                    + 3 * d1 * (4 + 3 * K) * (P034 * P041 + P033 * P042 - P031 * P044)
                    + 4 * (1 + 8 * K + K**2) * (P034 * P043 - P033 * P044)
                )
            )
            * S1**2
        )
        + B1**2
        * (
            C1**2
            * (
                2 * (9 + 12 * K + K**2) * P031 * P042
                - 8
                * d1**2
                * kj**3
                * (4 + 3 * K)
                * (P033 * P041 - P034 * P042 - P031 * P043)
                - 16 * d1**2 * kj**4 * (1 + 3 * K) * (P034 * P043 - P033 * P044)
                + kj
                * (9 + 12 * K + K**2)
                * (
                    P033 * (P041 - 3 * P042)
                    + P034 * (3 * P041 - P042)
                    - P031 * (P043 + 3 * P044)
                )
                + P032
                * (
                    2 * (-9 - 12 * K - K**2 + 2 * d1**2 * kj**2 * (7 + 3 * K)) * P041
                    + kj
                    * (
                        3 * (9 + 12 * K + K**2 + 2 * d1 * kj * (8 + 5 * K)) * P043
                        + (9 + 12 * K + K**2 - 8 * d1**2 * kj**2 * (4 + 3 * K)) * P044
                    )
                )
                - 2
                * kj**2
                * (
                    2 * d1**2 * (7 + 3 * K) * P031 * P042
                    + 3 * d1 * (8 + 5 * K) * (P034 * P041 + P033 * P042 - P031 * P044)
                    + 2 * (9 + 12 * K + K**2) * (P034 * P043 - P033 * P044)
                )
            )
            + 3
            * C1
            * (
                (-(36 + 21 * K + K**2)) * P031 * P042
                + 4 * kj**2 * K * (5 + K) * (P034 * P043 - P033 * P044)
                + P032
                * (
                    (36 + 21 * K + K**2) * P041
                    - 2 * kj * ((13 + 11 * K + K**2) * P043 + (3 + 2 * K) * P044)
                )
                + 2
                * kj
                * (
                    (-P034) * ((13 + 11 * K + K**2) * P041 - (3 + 2 * K) * P042)
                    + P033 * (-(3 + 2 * K) * P041 + (13 + 11 * K + K**2) * P042)
                    + P031 * ((3 + 2 * K) * P043 + (13 + 11 * K + K**2) * P044)
                )
            )
            * S1
            + (
                (106 + 35 * K + K**2) * P031 * P042
                + 8
                * d1**2
                * kj**3
                * (4 + 3 * K)
                * (P033 * P041 - P034 * P042 - P031 * P043)
                + 16 * d1**2 * kj**4 * (1 + 3 * K) * (P034 * P043 - P033 * P044)
                + kj
                * (
                    P034
                    * (3 * (9 + 12 * K + K**2) * P041 + (-17 + 2 * K + K**2) * P042)
                    - P033
                    * ((-17 + 2 * K + K**2) * P041 + 3 * (9 + 12 * K + K**2) * P042)
                    + P031
                    * ((-17 + 2 * K + K**2) * P043 - 3 * (9 + 12 * K + K**2) * P044)
                )
                - P032
                * (
                    (106 + 35 * K + K**2 + 4 * d1**2 * kj**2 * (7 + 3 * K)) * P041
                    + kj
                    * (
                        -3 * (9 + 12 * K + K**2 - 2 * d1 * kj * (8 + 5 * K)) * P043
                        + (-17 + 2 * K + K**2 - 8 * d1**2 * kj**2 * (4 + 3 * K)) * P044
                    )
                )
                + 2
                * kj**2
                * (
                    2 * d1**2 * (7 + 3 * K) * P031 * P042
                    + 3 * d1 * (8 + 5 * K) * (P034 * P041 + P033 * P042 - P031 * P044)
                    - 2 * (-1 + K + 2 * K**2) * (P034 * P043 - P033 * P044)
                )
            )
            * S1**2
        )
    )

    D7 = (
        -3
        * (C1 - S1)
        * (
            2
            * B1
            * (
                2
                * C1
                * (
                    2 * (3 + K) * P031 * P042
                    + 4 * d1**2 * kj**3 * (-P033 * P041 + P034 * P042 + P031 * P043)
                    + kj
                    * (3 + K)
                    * (
                        3 * P034 * P041
                        + P033 * (P041 - 3 * P042)
                        - P034 * P042
                        - P031 * P043
                        - 3 * P031 * P044
                    )
                    + 8 * d1**2 * kj**4 * (-P034 * P043 + P033 * P044)
                    + P032
                    * (
                        2 * (-3 + d1**2 * kj**2 - K) * P041
                        + kj
                        * (
                            3 * (3 + 2 * d1 * kj + K) * P043
                            + (3 - 4 * d1**2 * kj**2 + K) * P044
                        )
                    )
                    - 2
                    * kj**2
                    * (
                        d1**2 * P031 * P042
                        + 3 * d1 * (P034 * P041 + P033 * P042 - P031 * P044)
                        + 2 * (3 + K) * (P034 * P043 - P033 * P044)
                    )
                )
                + (
                    -3 * (8 + K) * P031 * P042
                    + 8 * d1**2 * kj**3 * (-P033 * P041 + P034 * P042 + P031 * P043)
                    - 16 * d1**2 * kj**4 * (P034 * P043 - P033 * P044)
                    + P032
                    * (
                        (4 * d1**2 * kj**2 + 3 * (8 + K)) * P041
                        - 2
                        * kj
                        * (
                            (9 - 6 * d1 * kj + 3 * K) * P043
                            + (3 + 4 * d1**2 * kj**2) * P044
                        )
                    )
                    - 6
                    * kj
                    * (
                        P034 * ((3 + K) * P041 - P042)
                        + P033 * (P041 - (3 + K) * P042)
                        - P031 * (P043 + (3 + K) * P044)
                    )
                    - 4
                    * kj**2
                    * (
                        d1**2 * P031 * P042
                        - 3 * K * P034 * P043
                        + 3 * K * P033 * P044
                        + 3 * d1 * (P034 * P041 + P033 * P042 - P031 * P044)
                    )
                )
                * S1
            )
            + B2
            * (
                C1
                * (
                    3 * (8 + K) * P031 * P042
                    + 16 * d1**2 * kj**3 * (P033 * P041 - P034 * P042 - P031 * P043)
                    + 32 * d1**2 * kj**4 * (P034 * P043 - P033 * P044)
                    + P032
                    * (
                        -(8 * d1**2 * kj**2 + 3 * (8 + K)) * P041
                        + 2
                        * kj
                        * (
                            3 * (3 - 4 * d1 * kj + K) * P043
                            + (3 + 8 * d1**2 * kj**2) * P044
                        )
                    )
                    + 6
                    * kj
                    * (
                        P034 * ((3 + K) * P041 - P042)
                        + P033 * (P041 - (3 + K) * P042)
                        - P031 * (P043 + (3 + K) * P044)
                    )
                    + 4
                    * kj**2
                    * (
                        2 * d1**2 * P031 * P042
                        - 3 * K * P034 * P043
                        + 3 * K * P033 * P044
                        + 6 * d1 * (P034 * P041 + P033 * P042 - P031 * P044)
                    )
                )
                + (
                    24
                    * d1
                    * kj**2
                    * (P034 * P041 + P033 * P042 - P032 * P043 - P031 * P044)
                    + 6
                    * kj
                    * (
                        P031 * P043
                        - 3 * P032 * P043
                        + P034 * (-3 * P041 + P042 + 8 * kj * P043)
                        + 3 * P031 * P044
                        - P032 * P044
                        - P033 * (P041 - 3 * P042 + 8 * kj * P044)
                    )
                    + K
                    * (
                        -5 * P031 * P042
                        + P032 * (5 * P041 - 6 * kj * P043 - 4 * kj * P044)
                        + kj
                        * (
                            -4 * P033 * P041
                            - 6 * P034 * P041
                            + 6 * P033 * P042
                            + 4 * P034 * P042
                            + 4 * P031 * P043
                            + 6 * P031 * P044
                        )
                        + 4 * kj**2 * (P034 * P043 - P033 * P044)
                    )
                    - 8
                    * d1**2
                    * kj**2
                    * (
                        -P031 * P042
                        + 2 * kj * (-P033 * P041 + P034 * P042 + P031 * P043)
                        + P032 * (P041 - 2 * kj * P044)
                        + kj**2 * (-4 * P034 * P043 + 4 * P033 * P044)
                    )
                )
                * S1
            )
        )
    )
    D8 = (
        9
        * (
            -2 * P031 * P042
            + kj
            * (
                -P033 * P041
                - 3 * P034 * P041
                + 3 * P033 * P042
                + P034 * P042
                + P031 * P043
                + 3 * P031 * P044
            )
            + 4 * kj**2 * (P034 * P043 - P033 * P044)
            + P032 * (2 * P041 - kj * (3 * P043 + P044))
        )
        * (C1 - S1) ** 2
    )

    return [D1, D2, D3, D4, D5, D6, D7, D8]


def get_num1_2visc(P0, Pzsf, d1, kj, B1, B2, K):

    P011, P012, P013, P014 = P0[0, :]
    P021, P022, P023, P024 = P0[1, :]
    P031, P032, P033, P034 = P0[2, :]
    P041, P042, P043, P044 = P0[3, :]

    # Calculating hyperbolic functions
    C1 = cosh(kj * d1)
    S1 = sinh(kj * d1)
    Pf1 = Pzsf[0]
    Pf2 = Pzsf[1]

    N1 = (
        64
        * B1**4
        * B2
        * K**3
        * (
            P012 * P041 * Pf1
            - P011 * P042 * Pf1
            - P012 * P031 * Pf2
            + P011 * P032 * Pf2
        )
        * (
            -2 * B1 * B2 * C1 * S1
            + B2**2 * (C1**2 * (1 + d1**2 * kj**2) - d1**2 * kj**2 * S1**2)
            + B1**2 * ((-(C1**2)) * d1**2 * kj**2 + (1 + d1**2 * kj**2) * S1**2)
        )
    )

    N2 = (
        -64
        * B1**3
        * K**2
        * (
            B1**3
            * (
                P012 * P041 * Pf1
                - P011 * P042 * Pf1
                - P012 * P031 * Pf2
                + P011 * P032 * Pf2
            )
            * (C1**2 * d1**2 * kj**2 - (1 + d1**2 * kj**2) * S1**2)
            + B1
            * B2**2
            * (
                C1**2
                * (
                    (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + d1**2 * kj**2 * (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        P014 * P041 * Pf1
                        - P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + P012
                    * (
                        (-(1 + d1**2 * kj**2)) * (2 + K) * P041 * Pf1
                        + (2 + K) * P031 * Pf2
                        + d1**2 * kj**2 * (2 + K) * P031 * Pf2
                        + kj * K * (P043 * Pf1 - P033 * Pf2)
                    )
                )
                + C1
                * (7 + 3 * K)
                * (
                    P012 * P041 * Pf1
                    - P011 * P042 * Pf1
                    - P012 * P031 * Pf2
                    + P011 * P032 * Pf2
                )
                * S1
                + kj
                * (
                    d1**2
                    * kj
                    * (2 + K)
                    * (
                        P012 * P041 * Pf1
                        - P011 * P042 * Pf1
                        - P012 * P031 * Pf2
                        + P011 * P032 * Pf2
                    )
                    + K
                    * (
                        P014 * P041 * Pf1
                        - P013 * P042 * Pf1
                        + P012 * P043 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        - P012 * P033 * Pf2
                        + P011 * P034 * Pf2
                    )
                )
                * S1**2
            )
            + B1**2
            * B2
            * (
                C1**2
                * d1
                * kj**2
                * (
                    K
                    * (
                        (-P014) * P041 * Pf1
                        - P013 * P042 * Pf1
                        + P012 * P043 * Pf1
                        + P011 * P044 * Pf1
                        + P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        - P012 * P033 * Pf2
                        - P011 * P034 * Pf2
                    )
                    + d1
                    * (
                        (-(3 + 2 * K)) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            (-P013) * P041 * Pf1
                            + P014 * P042 * Pf1
                            + P011 * P043 * Pf1
                            + P013 * P031 * Pf2
                            - P014 * P032 * Pf2
                            - P011 * P033 * Pf2
                        )
                        + P012
                        * (
                            (3 + 2 * K) * P041 * Pf1
                            - (3 + 2 * K) * P031 * Pf2
                            + kj * K * ((-P044) * Pf1 + P034 * Pf2)
                        )
                    )
                )
                + C1
                * (
                    (-(3 + K)) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        (-P014) * P041 * Pf1
                        + P013 * P042 * Pf1
                        + P011 * P044 * Pf1
                        + P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        - P011 * P034 * Pf2
                    )
                    + P012
                    * (
                        (3 + K) * P041 * Pf1
                        - (3 + K) * P031 * Pf2
                        + kj * K * ((-P043) * Pf1 + P033 * Pf2)
                    )
                )
                * S1
                + (
                    2 * (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + d1
                    * kj**2
                    * K
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (3 + 2 * K) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                    )
                    + P012
                    * (
                        (-(2 * (2 + K) + d1**2 * kj**2 * (3 + 2 * K))) * P041 * Pf1
                        + 2 * (2 + K) * P031 * Pf2
                        + d1 * kj**2 * K * ((-P043) * Pf1 + P033 * Pf2)
                        + d1**2
                        * kj**2
                        * (
                            (3 + 2 * K) * P031 * Pf2
                            + kj * K * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B2**3
            * (
                C1**2
                * (
                    (3 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + d1
                    * kj**2
                    * K
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                    )
                    + P012
                    * (
                        (-(3 + K + d1**2 * kj**2 * (2 + K))) * P041 * Pf1
                        + (3 + K) * P031 * Pf2
                        + d1 * kj**2 * K * ((-P043) * Pf1 + P033 * Pf2)
                        + d1**2
                        * kj**2
                        * ((2 + K) * P031 * Pf2 + kj * K * (P044 * Pf1 - P034 * Pf2))
                    )
                )
                + C1
                * kj
                * K
                * (
                    (-P014) * P041 * Pf1
                    + P013 * P042 * Pf1
                    - P012 * P043 * Pf1
                    + P011 * P044 * Pf1
                    + P014 * P031 * Pf2
                    - P013 * P032 * Pf2
                    + P012 * P033 * Pf2
                    - P011 * P034 * Pf2
                )
                * S1
                + d1
                * kj**2
                * (
                    K
                    * (
                        (-P014) * P041 * Pf1
                        - P013 * P042 * Pf1
                        + P012 * P043 * Pf1
                        + P011 * P044 * Pf1
                        + P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        - P012 * P033 * Pf2
                        - P011 * P034 * Pf2
                    )
                    + d1
                    * (
                        (-(2 + K)) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            (-P013) * P041 * Pf1
                            + P014 * P042 * Pf1
                            + P011 * P043 * Pf1
                            + P013 * P031 * Pf2
                            - P014 * P032 * Pf2
                            - P011 * P033 * Pf2
                        )
                        + P012
                        * (
                            (2 + K) * P041 * Pf1
                            - (2 + K) * P031 * Pf2
                            + kj * K * ((-P044) * Pf1 + P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
        )
    )

    N3 = (
        16
        * B1**2
        * K
        * (
            B1**2
            * B2
            * (
                C1**2
                * (
                    2
                    * d1
                    * kj**2
                    * K
                    * (7 + 3 * K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P012 * P043 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P012 * P033 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (12 + 14 * K + 5 * K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 6
                        * kj
                        * K
                        * (2 + K)
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        + 4
                        * kj**2
                        * K**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                        + P012
                        * (
                            (-(12 + 14 * K + 5 * K**2)) * P041 * Pf1
                            + (12 + 14 * K + 5 * K**2) * P031 * Pf2
                            + 6 * kj * K * (2 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                    + K
                    * (
                        (-(8 + K)) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 4
                        * kj**2
                        * K
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                        + P012
                        * (
                            (8 + K) * P041 * Pf1
                            - (8 + K) * P031 * Pf2
                            - 2
                            * kj
                            * (
                                (3 + K) * P043 * Pf1
                                + P044 * Pf1
                                - ((3 + K) * P033 + P034) * Pf2
                            )
                        )
                        - 2
                        * kj
                        * (
                            P014
                            * (
                                (3 + K) * P041 * Pf1
                                - P042 * Pf1
                                - 3 * P031 * Pf2
                                - K * P031 * Pf2
                                + P032 * Pf2
                            )
                            + P013
                            * (
                                P041 * Pf1
                                - (3 + K) * P042 * Pf1
                                + (-P031 + (3 + K) * P032) * Pf2
                            )
                            + P011
                            * (
                                (-P043) * Pf1
                                - (3 + K) * P044 * Pf1
                                + (P033 + (3 + K) * P034) * Pf2
                            )
                        )
                    )
                )
                - 2
                * C1
                * (7 + 3 * K)
                * (
                    (-(3 + K)) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        (-P014) * P041 * Pf1
                        + P013 * P042 * Pf1
                        + P011 * P044 * Pf1
                        + P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        - P011 * P034 * Pf2
                    )
                    + P012
                    * (
                        (3 + K) * P041 * Pf1
                        - (3 + K) * P031 * Pf2
                        + kj * K * ((-P043) * Pf1 + P033 * Pf2)
                    )
                )
                * S1
                + (
                    -2
                    * d1
                    * kj**2
                    * K
                    * (7 + 3 * K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    - 2
                    * (
                        (9 + 16 * K + 3 * K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (3 + K)
                        * (
                            P014 * P041 * Pf1
                            - P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            + P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                    )
                    + d1**2
                    * kj**2
                    * (
                        (-(12 + 14 * K + 5 * K**2)) * P011 * (P042 * Pf1 - P032 * Pf2)
                        - 6
                        * kj
                        * K
                        * (2 + K)
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        - 4
                        * kj**2
                        * K**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        (
                            18
                            + 32 * K
                            + 6 * K**2
                            + d1**2 * kj**2 * (12 + 14 * K + 5 * K**2)
                        )
                        * P041
                        * Pf1
                        - 2 * (9 + 16 * K + 3 * K**2) * P031 * Pf2
                        - 2 * kj * K * (3 + K) * (P043 * Pf1 - P033 * Pf2)
                        - 6 * d1**2 * kj**3 * K * (2 + K) * (P044 * Pf1 - P034 * Pf2)
                        - d1
                        * kj**2
                        * (
                            12 * d1 * P031 * Pf2
                            + K**2
                            * (-6 * P043 * Pf1 + 5 * d1 * P031 * Pf2 + 6 * P033 * Pf2)
                            - 14 * K * (P043 * Pf1 - (d1 * P031 + P033) * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + 2
            * B1
            * B2**2
            * (
                2
                * C1**2
                * (
                    (-d1)
                    * kj**2
                    * K
                    * (2 + K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    - d1**2
                    * kj**2
                    * (2 + K)
                    * (
                        (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * ( 
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                    )
                    - (3 + K)
                    * (
                        (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P014 * P041 * Pf1
                            - P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            + P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        (2 + K) * (3 + K + d1**2 * kj**2 * (2 + K)) * P041 * Pf1
                        - (6 + 5 * K + K**2) * P031 * Pf2
                        - kj * K * (3 + K) * (P043 * Pf1 - P033 * Pf2)
                        - d1**2 * kj**3 * K * (2 + K) * (P044 * Pf1 - P034 * Pf2)
                        - d1
                        * kj**2
                        * (2 + K)
                        * (
                            2 * d1 * P031 * Pf2
                            + K * ((-P043) * Pf1 + (d1 * P031 + P033) * Pf2)
                        )
                    )
                )
                + C1
                * (
                    3 * (4 + 7 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * kj**2
                    * K**2
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + P012
                    * (
                        -3 * (4 + 7 * K + K**2) * P041 * Pf1
                        + 3 * (4 + 7 * K + K**2) * P031 * Pf2
                        + 2
                        * kj
                        * K
                        * (
                            (2 + K) * P043 * Pf1
                            + P044 * Pf1
                            - ((2 + K) * P033 + P034) * Pf2
                        )
                    )
                    - 2
                    * kj
                    * K
                    * (
                        P014
                        * (
                            (-(2 + K)) * P041 * Pf1
                            + P042 * Pf1
                            + ((2 + K) * P031 - P032) * Pf2
                        )
                        + P013
                        * (
                            (-P041) * Pf1
                            + (2 + K) * P042 * Pf1
                            + (P031 - (2 + K) * P032) * Pf2
                        )
                        + P011
                        * (
                            P043 * Pf1
                            + (2 + K) * P044 * Pf1
                            - (P033 + (2 + K) * P034) * Pf2
                        )
                    )
                )
                * S1
                + 2
                * kj
                * (
                    (-K)
                    * (3 + K)
                    * (
                        P014 * P041 * Pf1
                        - P013 * P042 * Pf1
                        + P012 * P043 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        - P012 * P033 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1
                    * kj
                    * K
                    * (2 + K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P012 * P043 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P012 * P033 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1**2
                    * kj
                    * (2 + K)
                    * (
                        (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        + P012
                        * (
                            (-(2 + K)) * P041 * Pf1
                            + (2 + K) * P031 * Pf2
                            + kj * K * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + 2
            * B1**3
            * (
                C1**2
                * d1
                * kj**2
                * (
                    2
                    * K
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P012 * P043 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P012 * P033 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1
                    * (
                        (4 + 5 * K) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 2
                        * kj
                        * K
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        + P012
                        * (
                            (-(4 + 5 * K)) * P041 * Pf1
                            + (4 + 5 * K) * P031 * Pf2
                            + 2 * kj * K * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                + C1
                * K
                * (
                    3 * P011 * (P042 * Pf1 - P032 * Pf2)
                    + P012
                    * (
                        -3 * P041 * Pf1
                        + 2 * kj * P043 * Pf1
                        + 3 * P031 * Pf2
                        - 2 * kj * P033 * Pf2
                    )
                    + 2
                    * kj
                    * (
                        P014 * P041 * Pf1
                        - P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                )
                * S1
                + (
                    -4 * (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * d1
                    * kj**2
                    * K
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (-(4 + 5 * K)) * P011 * (P042 * Pf1 - P032 * Pf2)
                        - 2
                        * kj
                        * K
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                    )
                    + P012
                    * (
                        (4 * (2 + K) + d1**2 * kj**2 * (4 + 5 * K)) * P041 * Pf1
                        - 4 * (2 + K) * P031 * Pf2
                        + 2 * d1 * kj**2 * K * (P043 * Pf1 - P033 * Pf2)
                        - d1**2
                        * kj**2
                        * (
                            (4 + 5 * K) * P031 * Pf2
                            + 2 * kj * K * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B2**3
            * (
                C1**2
                * (
                    (-(9 + 12 * K + K**2)) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * d1
                    * kj**2
                    * K
                    * (5 + K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (-(4 + 8 * K + K**2)) * P011 * (P042 * Pf1 - P032 * Pf2)
                        - 2
                        * kj
                        * K
                        * (4 + K)
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        - 4
                        * kj**2
                        * K**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        (9 + 12 * K + K**2 + d1**2 * kj**2 * (4 + 8 * K + K**2))
                        * P041
                        * Pf1
                        - (9 + 12 * K + K**2) * P031 * Pf2
                        + 2 * d1 * kj**2 * K * (5 + K) * (P043 * Pf1 - P033 * Pf2)
                        - d1**2
                        * kj**2
                        * (
                            (4 + 8 * K + K**2) * P031 * Pf2
                            + 2 * kj * K * (4 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                + 2
                * C1
                * kj
                * K
                * (5 + K)
                * (
                    P014 * P041 * Pf1
                    - P013 * P042 * Pf1
                    + P012 * P043 * Pf1
                    - P011 * P044 * Pf1
                    - P014 * P031 * Pf2
                    + P013 * P032 * Pf2
                    - P012 * P033 * Pf2
                    + P011 * P034 * Pf2
                )
                * S1
                + (
                    P011 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * kj
                    * K
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    + 2
                    * d1**2
                    * kj**3
                    * K
                    * (4 + K)
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    + 4
                    * d1**2
                    * kj**4
                    * K**2
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + kj**2
                    * (
                        d1**2 * (4 + 8 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 2
                        * d1
                        * K
                        * (5 + K)
                        * (
                            P014 * P041 * Pf1
                            + P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            - P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                        + 4
                        * K**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        (-(1 + d1**2 * kj**2 * (4 + 8 * K + K**2))) * P041 * Pf1
                        - 2 * kj * K * P044 * Pf1
                        + P031 * Pf2
                        + 2 * kj * K * P034 * Pf2
                        - 2 * d1 * kj**2 * K * (5 + K) * (P043 * Pf1 - P033 * Pf2)
                        + d1**2
                        * kj**2
                        * (
                            (4 + 8 * K + K**2) * P031 * Pf2
                            + 2 * kj * K * (4 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
        )
    )

    N4 = (
        16
        * B1
        * (
            B1
            * B2**2
            * (
                C1**2
                * (
                    -2
                    * d1
                    * kj**2
                    * K
                    * (10 + 7 * K + K**2)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    - (9 + 12 * K + K**2)
                    * (
                        (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P014 * P041 * Pf1
                            - P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            + P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                    )
                    + d1**2
                    * kj**2
                    * (2 + K)
                    * (
                        (-(4 + 8 * K + K**2)) * P011 * (P042 * Pf1 - P032 * Pf2)
                        - 2
                        * kj
                        * K
                        * (4 + K)
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        - 4
                        * kj**2
                        * K**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        (2 + K)
                        * (9 + 12 * K + K**2 + d1**2 * kj**2 * (4 + 8 * K + K**2))
                        * P041
                        * Pf1
                        - (18 + 33 * K + 14 * K**2 + K**3) * P031 * Pf2
                        - kj * K * (9 + 12 * K + K**2) * (P043 * Pf1 - P033 * Pf2)
                        - 2
                        * d1**2
                        * kj**3
                        * K
                        * (8 + 6 * K + K**2)
                        * (P044 * Pf1 - P034 * Pf2)
                        - d1
                        * kj**2
                        * (2 + K)
                        * (
                            4 * d1 * P031 * Pf2
                            + K**2
                            * (-2 * P043 * Pf1 + d1 * P031 * Pf2 + 2 * P033 * Pf2)
                            + 2
                            * K
                            * (-5 * P043 * Pf1 + 4 * d1 * P031 * Pf2 + 5 * P033 * Pf2)
                        )
                    )
                )
                + C1
                * K
                * (
                    (36 + 21 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * kj**2
                    * K
                    * (5 + K)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + P012
                    * (
                        (-(36 + 21 * K + K**2)) * P041 * Pf1
                        + (36 + 21 * K + K**2) * P031 * Pf2
                        + 2
                        * kj
                        * (
                            (10 + 7 * K + K**2) * P043 * Pf1
                            + (3 + 2 * K) * P044 * Pf1
                            - ((10 + 7 * K + K**2) * P033 + (3 + 2 * K) * P034) * Pf2
                        )
                    )
                    + 2
                    * kj
                    * (
                        P014
                        * (
                            (10 + 7 * K + K**2) * P041 * Pf1
                            - (3 + 2 * K) * P042 * Pf1
                            - ((10 + 7 * K + K**2) * P031 - (3 + 2 * K) * P032) * Pf2
                        )
                        + P013
                        * (
                            (3 + 2 * K) * P041 * Pf1
                            - (10 + 7 * K + K**2) * P042 * Pf1
                            + ((-(3 + 2 * K)) * P031 + (10 + 7 * K + K**2) * P032) * Pf2
                        )
                        + P011
                        * (
                            (-(3 + 2 * K)) * P043 * Pf1
                            - (10 + 7 * K + K**2) * P044 * Pf1
                            + ((3 + 2 * K) * P033 + (10 + 7 * K + K**2) * P034) * Pf2
                        )
                    )
                )
                * S1
                + (
                    (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 2
                    * d1**2
                    * kj**3
                    * K
                    * (8 + 6 * K + K**2)
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    + 4
                    * d1**2
                    * kj**4
                    * K**2
                    * (2 + K)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + kj**2
                    * (2 + K)
                    * (
                        d1**2 * (4 + 8 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 2
                        * d1
                        * K
                        * (5 + K)
                        * (
                            P014 * P041 * Pf1
                            + P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            - P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                        + 4
                        * K**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                    )
                    + kj
                    * K
                    * (
                        P014
                        * (
                            (-(9 + 12 * K + K**2)) * P041 * Pf1
                            + 2 * (2 + K) * P042 * Pf1
                            + ((9 + 12 * K + K**2) * P031 - 2 * (2 + K) * P032) * Pf2
                        )
                        + P013
                        * (
                            -2 * (2 + K) * P041 * Pf1
                            + (9 + 12 * K + K**2) * P042 * Pf1
                            + (2 * (2 + K) * P031 - (9 + 12 * K + K**2) * P032) * Pf2
                        )
                        + P011
                        * (
                            2 * (2 + K) * P043 * Pf1
                            + (9 + 12 * K + K**2) * P044 * Pf1
                            - (2 * (2 + K) * P033 + (9 + 12 * K + K**2) * P034) * Pf2
                        )
                    )
                    + P012
                    * (
                        (-(2 + K))
                        * (1 + d1**2 * kj**2 * (4 + 8 * K + K**2))
                        * P041
                        * Pf1
                        + (2 + K) * P031 * Pf2
                        + 2
                        * d1**2
                        * kj**3
                        * K
                        * (8 + 6 * K + K**2)
                        * (P044 * Pf1 - P034 * Pf2)
                        + kj
                        * K
                        * (
                            (-(9 + 12 * K + K**2)) * P043 * Pf1
                            - 2 * (2 + K) * P044 * Pf1
                            + ((9 + 12 * K + K**2) * P033 + 2 * (2 + K) * P034) * Pf2
                        )
                        + d1
                        * kj**2
                        * (2 + K)
                        * (
                            4 * d1 * P031 * Pf2
                            + K**2
                            * (-2 * P043 * Pf1 + d1 * P031 * Pf2 + 2 * P033 * Pf2)
                            + 2
                            * K
                            * (-5 * P043 * Pf1 + 4 * d1 * P031 * Pf2 + 5 * P033 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B1**2
            * B2
            * (
                C1**2
                * (
                    d1
                    * kj**2
                    * K
                    * (16 + 11 * K + 2 * K**2)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P012 * P043 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P012 * P033 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (4 + 4 * K + 5 * K**2 + K**3) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 2
                        * kj
                        * K
                        * (6 + 4 * K + K**2)
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        + 4
                        * kj**2
                        * K**2
                        * (3 + K)
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                        + P012
                        * (
                            (-(4 + 4 * K + 5 * K**2 + K**3)) * P041 * Pf1
                            + (4 + 4 * K + 5 * K**2 + K**3) * P031 * Pf2
                            + 2
                            * kj
                            * K
                            * (6 + 4 * K + K**2)
                            * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                    + K
                    * (3 + K)
                    * (
                        (-(8 + K)) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 4
                        * kj**2
                        * K
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                        + P012
                        * (
                            (8 + K) * P041 * Pf1
                            - (8 + K) * P031 * Pf2
                            - 2
                            * kj
                            * (
                                (3 + K) * P043 * Pf1
                                + P044 * Pf1
                                - ((3 + K) * P033 + P034) * Pf2
                            )
                        )
                        - 2
                        * kj
                        * (
                            P014
                            * (
                                (3 + K) * P041 * Pf1
                                - P042 * Pf1
                                - 3 * P031 * Pf2
                                - K * P031 * Pf2
                                + P032 * Pf2
                            )
                            + P013
                            * (
                                P041 * Pf1
                                - (3 + K) * P042 * Pf1
                                + (-P031 + (3 + K) * P032) * Pf2
                            )
                            + P011
                            * (
                                (-P043) * Pf1
                                - (3 + K) * P044 * Pf1
                                + (P033 + (3 + K) * P034) * Pf2
                            )
                        )
                    )
                )
                + C1
                * (
                    3
                    * (12 + 25 * K + 10 * K**2 + K**3)
                    * P011
                    * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * kj**2
                    * K**2
                    * (3 + K)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + P012
                    * (
                        -3 * (12 + 25 * K + 10 * K**2 + K**3) * P041 * Pf1
                        + 3 * (12 + 25 * K + 10 * K**2 + K**3) * P031 * Pf2
                        + kj
                        * K
                        * (
                            (14 + 29 * K + 4 * K**2) * P043 * Pf1
                            + 2 * (3 + K) * P044 * Pf1
                            - ((14 + 29 * K + 4 * K**2) * P033 + 2 * (3 + K) * P034)
                            * Pf2
                        )
                    )
                    + kj
                    * K
                    * (
                        P014
                        * (
                            (14 + 29 * K + 4 * K**2) * P041 * Pf1
                            - 2 * (3 + K) * P042 * Pf1
                            - ((14 + 29 * K + 4 * K**2) * P031 - 2 * (3 + K) * P032)
                            * Pf2
                        )
                        + P013
                        * (
                            2 * (3 + K) * P041 * Pf1
                            - (14 + 29 * K + 4 * K**2) * P042 * Pf1
                            + (-2 * (3 + K) * P031 + (14 + 29 * K + 4 * K**2) * P032)
                            * Pf2
                        )
                        + P011
                        * (
                            -2 * (3 + K) * P043 * Pf1
                            - (14 + 29 * K + 4 * K**2) * P044 * Pf1
                            + (2 * (3 + K) * P033 + (14 + 29 * K + 4 * K**2) * P034)
                            * Pf2
                        )
                    )
                )
                * S1
                + (
                    (-(4 + 35 * K + 24 * K**2 + 2 * K**3))
                    * P011
                    * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * d1**2
                    * kj**3
                    * K
                    * (6 + 4 * K + K**2)
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    - 4
                    * d1**2
                    * kj**4
                    * K**2
                    * (3 + K)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + d1
                    * kj**2
                    * (
                        4 * d1 * P011 * ((-P042) * Pf1 + P032 * Pf2)
                        + K**2
                        * (
                            -11 * P014 * P041 * Pf1
                            - 5 * d1 * P011 * P042 * Pf1
                            - 11 * P013 * P042 * Pf1
                            + 11 * P011 * P044 * Pf1
                            + 11 * P014 * P031 * Pf2
                            + 5 * d1 * P011 * P032 * Pf2
                            + 11 * P013 * P032 * Pf2
                            - 11 * P011 * P034 * Pf2
                        )
                        + K**3
                        * (
                            -2 * P014 * P041 * Pf1
                            - d1 * P011 * P042 * Pf1
                            - 2 * P013 * P042 * Pf1
                            + 2 * P011 * P044 * Pf1
                            + 2 * P014 * P031 * Pf2
                            + d1 * P011 * P032 * Pf2
                            + 2 * P013 * P032 * Pf2
                            - 2 * P011 * P034 * Pf2
                        )
                        - 4
                        * K
                        * (
                            4 * P014 * P041 * Pf1
                            + d1 * P011 * P042 * Pf1
                            + 4 * P013 * P042 * Pf1
                            - 4 * P011 * P044 * Pf1
                            - 4 * P014 * P031 * Pf2
                            - d1 * P011 * P032 * Pf2
                            - 4 * P013 * P032 * Pf2
                            + 4 * P011 * P034 * Pf2
                        )
                    )
                    + kj
                    * K
                    * (
                        P014
                        * (
                            -2 * (3 + K) ** 2 * P041 * Pf1
                            + (-2 + K) * P042 * Pf1
                            + (2 * (3 + K) ** 2 * P031 - (-2 + K) * P032) * Pf2
                        )
                        + P013
                        * (
                            (-(-2 + K)) * P041 * Pf1
                            + 2 * (3 + K) ** 2 * P042 * Pf1
                            + ((-2 + K) * P031 - 2 * (3 + K) ** 2 * P032) * Pf2
                        )
                        + P011
                        * (
                            (-2 + K) * P043 * Pf1
                            + 2 * (3 + K) ** 2 * P044 * Pf1
                            - ((-2 + K) * P033 + 2 * (3 + K) ** 2 * P034) * Pf2
                        )
                    )
                    + P012
                    * (
                        (
                            4
                            + 35 * K
                            + 24 * K**2
                            + 2 * K**3
                            + d1**2 * kj**2 * (4 + 4 * K + 5 * K**2 + K**3)
                        )
                        * P041
                        * Pf1
                        - (4 + 35 * K + 24 * K**2 + 2 * K**3) * P031 * Pf2
                        - 2
                        * d1**2
                        * kj**3
                        * K
                        * (6 + 4 * K + K**2)
                        * (P044 * Pf1 - P034 * Pf2)
                        - kj
                        * K
                        * (
                            2 * (3 + K) ** 2 * P043 * Pf1
                            + (-2 + K) * P044 * Pf1
                            - (2 * (3 + K) ** 2 * P033 + (-2 + K) * P034) * Pf2
                        )
                        - d1
                        * kj**2
                        * (
                            4 * d1 * P031 * Pf2
                            + K**3
                            * (-2 * P043 * Pf1 + d1 * P031 * Pf2 + 2 * P033 * Pf2)
                            + 4
                            * K
                            * (-4 * P043 * Pf1 + d1 * P031 * Pf2 + 4 * P033 * Pf2)
                            + K**2
                            * (-11 * P043 * Pf1 + 5 * d1 * P031 * Pf2 + 11 * P033 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B2**3
            * K
            * (
                C1**2
                * (
                    -3 * (3 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - d1
                    * kj**2
                    * (6 + 5 * K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + 2
                    * d1**2
                    * kj**2
                    * (
                        (-(2 + K)) * P011 * (P042 * Pf1 - P032 * Pf2)
                        - 2
                        * kj
                        * (1 + K)
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        - 4
                        * kj**2
                        * K
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        (2 * d1**2 * kj**2 * (2 + K) + 3 * (3 + K)) * P041 * Pf1
                        - 3 * (3 + K) * P031 * Pf2
                        + d1 * kj**2 * (6 + 5 * K) * (P043 * Pf1 - P033 * Pf2)
                        - 2
                        * d1**2
                        * kj**2
                        * (
                            (2 + K) * P031 * Pf2
                            + 2 * kj * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                + C1
                * kj
                * (6 + 5 * K)
                * (
                    P014 * P041 * Pf1
                    - P013 * P042 * Pf1
                    + P012 * P043 * Pf1
                    - P011 * P044 * Pf1
                    - P014 * P031 * Pf2
                    + P013 * P032 * Pf2
                    - P012 * P033 * Pf2
                    + P011 * P034 * Pf2
                )
                * S1
                + (
                    P011 * (P042 * Pf1 - P032 * Pf2)
                    + 4
                    * d1**2
                    * kj**3
                    * (1 + K)
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    - kj
                    * (2 + K)
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    + 8
                    * d1**2
                    * kj**4
                    * K
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + kj**2
                    * (
                        2 * d1**2 * (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + d1
                        * (6 + 5 * K)
                        * (
                            P014 * P041 * Pf1
                            + P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            - P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                        + 8
                        * K
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        (-(1 + 2 * d1**2 * kj**2 * (2 + K))) * P041 * Pf1
                        - 2 * kj * P044 * Pf1
                        - kj * K * P044 * Pf1
                        + P031 * Pf2
                        + 2 * kj * P034 * Pf2
                        + kj * K * P034 * Pf2
                        - d1 * kj**2 * (6 + 5 * K) * (P043 * Pf1 - P033 * Pf2)
                        + 2
                        * d1**2
                        * kj**2
                        * (
                            (2 + K) * P031 * Pf2
                            + 2 * kj * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B1**3
            * (
                C1**2
                * (
                    2
                    * d1
                    * kj**2
                    * K
                    * (5 + 4 * K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P012 * P043 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P012 * P033 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (4 + 20 * K + 7 * K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 8
                        * kj
                        * K
                        * (1 + K)
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        + 4
                        * kj**2
                        * K**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                        + P012
                        * (
                            (-(4 + 20 * K + 7 * K**2)) * P041 * Pf1
                            + (4 + 20 * K + 7 * K**2) * P031 * Pf2
                            + 8 * kj * K * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                    + K**2
                    * (
                        -2 * P011 * P042 * Pf1
                        + 2 * P011 * P032 * Pf2
                        + 4
                        * kj**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                        + P012
                        * (
                            2 * P041 * Pf1
                            - 2 * P031 * Pf2
                            + kj
                            * (
                                -3 * P043 * Pf1
                                - P044 * Pf1
                                + 3 * P033 * Pf2
                                + P034 * Pf2
                            )
                        )
                        + kj
                        * (
                            P013
                            * ((-P041) * Pf1 + 3 * P042 * Pf1 + (P031 - 3 * P032) * Pf2)
                            + P014
                            * (
                                -3 * P041 * Pf1
                                + P042 * Pf1
                                + 3 * P031 * Pf2
                                - P032 * Pf2
                            )
                            + P011
                            * (P043 * Pf1 + 3 * P044 * Pf1 - (P033 + 3 * P034) * Pf2)
                        )
                    )
                )
                + C1
                * K
                * (7 + 3 * K)
                * (
                    3 * P011 * (P042 * Pf1 - P032 * Pf2)
                    + P012
                    * (
                        -3 * P041 * Pf1
                        + 2 * kj * P043 * Pf1
                        + 3 * P031 * Pf2
                        - 2 * kj * P033 * Pf2
                    )
                    + 2
                    * kj
                    * (
                        P014 * P041 * Pf1
                        - P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                )
                * S1
                + (
                    (-(16 + 33 * K + 6 * K**2)) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 8
                    * d1**2
                    * kj**3
                    * K
                    * (1 + K)
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    - 4
                    * d1**2
                    * kj**4
                    * K**2
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + d1
                    * kj**2
                    * (
                        4 * d1 * P011 * ((-P042) * Pf1 + P032 * Pf2)
                        + K**2
                        * (
                            -8 * P014 * P041 * Pf1
                            - 7 * d1 * P011 * P042 * Pf1
                            - 8 * P013 * P042 * Pf1
                            + 8 * P011 * P044 * Pf1
                            + 8 * P014 * P031 * Pf2
                            + 7 * d1 * P011 * P032 * Pf2
                            + 8 * P013 * P032 * Pf2
                            - 8 * P011 * P034 * Pf2
                        )
                        - 10
                        * K
                        * (
                            P014 * P041 * Pf1
                            + 2 * d1 * P011 * P042 * Pf1
                            + P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            - 2 * d1 * P011 * P032 * Pf2
                            - P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                    )
                    + kj
                    * K
                    * (
                        2
                        * (
                            P014 * P042 * Pf1
                            + P011 * P043 * Pf1
                            - P014 * P032 * Pf2
                            - P011 * P033 * Pf2
                        )
                        + P013
                        * (
                            (-2 + K) * P041 * Pf1
                            + 2 * P031 * Pf2
                            + K * (3 * P042 * Pf1 - (P031 + 3 * P032) * Pf2)
                        )
                        + K
                        * (
                            P014
                            * (
                                -3 * P041 * Pf1
                                - P042 * Pf1
                                + 3 * P031 * Pf2
                                + P032 * Pf2
                            )
                            + P011
                            * (
                                (-P043) * Pf1
                                + 3 * P044 * Pf1
                                + P033 * Pf2
                                - 3 * P034 * Pf2
                            )
                        )
                    )
                    + P012
                    * (
                        (
                            16
                            + 33 * K
                            + 6 * K**2
                            + d1**2 * kj**2 * (4 + 20 * K + 7 * K**2)
                        )
                        * P041
                        * Pf1
                        - 3 * kj * K**2 * P043 * Pf1
                        - 2 * kj * K * P044 * Pf1
                        + kj * K**2 * P044 * Pf1
                        - 16 * P031 * Pf2
                        - 33 * K * P031 * Pf2
                        - 6 * K**2 * P031 * Pf2
                        + 3 * kj * K**2 * P033 * Pf2
                        + 2 * kj * K * P034 * Pf2
                        - kj * K**2 * P034 * Pf2
                        + 2 * d1 * kj**2 * K * (5 + 4 * K) * (P043 * Pf1 - P033 * Pf2)
                        - d1**2
                        * kj**2
                        * (
                            (4 + 20 * K + 7 * K**2) * P031 * Pf2
                            + 8 * kj * K * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
        )
    )

    N5 = 4 * (
        B2**3
        * K
        * (
            C1**2
            * (
                -9 * P011 * P042 * Pf1
                + 9 * P011 * P032 * Pf2
                - 12
                * d1
                * kj**2
                * (
                    P014 * P041 * Pf1
                    + P013 * P042 * Pf1
                    - P011 * P044 * Pf1
                    - P014 * P031 * Pf2
                    - P013 * P032 * Pf2
                    + P011 * P034 * Pf2
                )
                + P012
                * (
                    (9 + 4 * d1**2 * kj**2) * P041 * Pf1
                    - 9 * P031 * Pf2
                    + 12 * d1 * kj**2 * (P043 * Pf1 - P033 * Pf2)
                    - 4
                    * d1**2
                    * kj**2
                    * (2 * kj * P044 * Pf1 + P031 * Pf2 - 2 * kj * P034 * Pf2)
                )
                + 4
                * d1**2
                * kj**2
                * (
                    (-P011) * P042 * Pf1
                    + P011 * P032 * Pf2
                    - 2
                    * kj
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    - 4
                    * kj**2
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
            )
            + 12
            * C1
            * kj
            * (
                P014 * P041 * Pf1
                - P013 * P042 * Pf1
                + P012 * P043 * Pf1
                - P011 * P044 * Pf1
                - P014 * P031 * Pf2
                + P013 * P032 * Pf2
                - P012 * P033 * Pf2
                + P011 * P034 * Pf2
            )
            * S1
            + (
                P011 * (P042 * Pf1 - P032 * Pf2)
                - 4
                * kj
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                + 8
                * d1**2
                * kj**3
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                + 16
                * d1**2
                * kj**4
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + P012
                * (
                    -4 * kj * P044 * Pf1
                    - P041 * (Pf1 + 4 * d1**2 * kj**2 * Pf1)
                    + P031 * Pf2
                    + 4 * kj * P034 * Pf2
                    - 12 * d1 * kj**2 * (P043 * Pf1 - P033 * Pf2)
                    + 4
                    * d1**2
                    * kj**2
                    * (2 * kj * P044 * Pf1 + P031 * Pf2 - 2 * kj * P034 * Pf2)
                )
                + 4
                * kj**2
                * (
                    d1**2 * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + 4
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
            )
            * S1**2
        )
        + B1**2
        * B2
        * (
            C1**2
            * (
                (-(72 + 105 * K + 20 * K**2 + K**3)) * P011 * (P042 * Pf1 - P032 * Pf2)
                + 4
                * d1**2
                * kj**3
                * (4 - 2 * K + K**2)
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                + 8
                * d1**2
                * kj**4
                * K
                * (6 + K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + 2
                * kj**2
                * (
                    d1**2 * (-4 - 10 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 2
                    * d1
                    * (6 - K + 2 * K**2)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + 2
                    * K
                    * (9 + 12 * K + K**2)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + P012
                * (
                    (
                        72
                        + 105 * K
                        + 20 * K**2
                        + K**3
                        - 2 * d1**2 * kj**2 * (-4 - 10 * K + K**2)
                    )
                    * P041
                    * Pf1
                    - (72 + 105 * K + 20 * K**2 + K**3) * P031 * Pf2
                    + 2
                    * d1
                    * kj**2
                    * (
                        -2 * (6 - K + 2 * K**2) * P043 * Pf1
                        + (
                            d1 * (-4 - 10 * K + K**2) * P031
                            + 2 * (6 - K + 2 * K**2) * P033
                        )
                        * Pf2
                    )
                    + 4 * d1**2 * kj**3 * (4 - 2 * K + K**2) * (P044 * Pf1 - P034 * Pf2)
                    - 2
                    * kj
                    * (9 + 12 * K + K**2)
                    * (
                        (3 + K) * P043 * Pf1
                        + P044 * Pf1
                        - ((3 + K) * P033 + P034) * Pf2
                    )
                )
                + 2
                * kj
                * (9 + 12 * K + K**2)
                * (
                    P014
                    * (
                        (-(3 + K)) * P041 * Pf1
                        + P042 * Pf1
                        + ((3 + K) * P031 - P032) * Pf2
                    )
                    + P013
                    * (
                        (-P041) * Pf1
                        + (3 + K) * P042 * Pf1
                        + (P031 - (3 + K) * P032) * Pf2
                    )
                    + P011
                    * (
                        P043 * Pf1
                        + (3 + K) * P044 * Pf1
                        - (P033 + (3 + K) * P034) * Pf2
                    )
                )
            )
            - 2
            * C1
            * (
                (-(108 + 99 * K + 24 * K**2 + K**3)) * P011 * (P042 * Pf1 - P032 * Pf2)
                + 4
                * kj**2
                * K
                * (15 + 8 * K + K**2)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + P012
                * (
                    (108 + 99 * K + 24 * K**2 + K**3) * P041 * Pf1
                    - (108 + 99 * K + 24 * K**2 + K**3) * P031 * Pf2
                    - 2
                    * kj
                    * (
                        (3 + 40 * K + 17 * K**2 + K**3) * P043 * Pf1
                        + (9 + 9 * K + 2 * K**2) * P044 * Pf1
                        - (
                            (3 + 40 * K + 17 * K**2 + K**3) * P033
                            + (9 + 9 * K + 2 * K**2) * P034
                        )
                        * Pf2
                    )
                )
                - 2
                * kj
                * (
                    P014
                    * (
                        (3 + 40 * K + 17 * K**2 + K**3) * P041 * Pf1
                        - (9 + 9 * K + 2 * K**2) * P042 * Pf1
                        - (
                            (3 + 40 * K + 17 * K**2 + K**3) * P031
                            - (9 + 9 * K + 2 * K**2) * P032
                        )
                        * Pf2
                    )
                    + P013
                    * (
                        (9 + 9 * K + 2 * K**2) * P041 * Pf1
                        - (3 + 40 * K + 17 * K**2 + K**3) * P042 * Pf1
                        - (
                            (9 + 9 * K + 2 * K**2) * P031
                            - (3 + 40 * K + 17 * K**2 + K**3) * P032
                        )
                        * Pf2
                    )
                    + P011
                    * (
                        (-(9 + 9 * K + 2 * K**2)) * P043 * Pf1
                        - (3 + 40 * K + 17 * K**2 + K**3) * P044 * Pf1
                        + (
                            (9 + 9 * K + 2 * K**2) * P033
                            + (3 + 40 * K + 17 * K**2 + K**3) * P034
                        )
                        * Pf2
                    )
                )
            )
            * S1
            + (
                (-(16 + 101 * K + 32 * K**2 + K**3)) * P011 * (P042 * Pf1 - P032 * Pf2)
                - 4
                * d1**2
                * kj**3
                * (4 - 2 * K + K**2)
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                - 8
                * d1**2
                * kj**4
                * K
                * (6 + K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + kj**2
                * (
                    -2 * d1**2 * (-4 - 10 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * d1
                    * (6 - K + 2 * K**2)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + 4
                    * K
                    * (1 + 8 * K + K**2)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + P012
                * (
                    (
                        16
                        + 101 * K
                        + 32 * K**2
                        + K**3
                        + 2 * d1**2 * kj**2 * (-4 - 10 * K + K**2)
                    )
                    * P041
                    * Pf1
                    - (16 + 101 * K + 32 * K**2 + K**3) * P031 * Pf2
                    - 2
                    * d1
                    * kj**2
                    * (
                        -2 * (6 - K + 2 * K**2) * P043 * Pf1
                        + (
                            d1 * (-4 - 10 * K + K**2) * P031
                            + 2 * (6 - K + 2 * K**2) * P033
                        )
                        * Pf2
                    )
                    - 4 * d1**2 * kj**3 * (4 - 2 * K + K**2) * (P044 * Pf1 - P034 * Pf2)
                    - 2
                    * kj
                    * (
                        (27 + 45 * K + 15 * K**2 + K**3) * P043 * Pf1
                        + (-7 + 8 * K + 3 * K**2) * P044 * Pf1
                        - (
                            (27 + 45 * K + 15 * K**2 + K**3) * P033
                            + (-7 + 8 * K + 3 * K**2) * P034
                        )
                        * Pf2
                    )
                )
                - 2
                * kj
                * (
                    P014
                    * (
                        (27 + 45 * K + 15 * K**2 + K**3) * P041 * Pf1
                        - (-7 + 8 * K + 3 * K**2) * P042 * Pf1
                        - (
                            (27 + 45 * K + 15 * K**2 + K**3) * P031
                            + (7 - 8 * K - 3 * K**2) * P032
                        )
                        * Pf2
                    )
                    + P013
                    * (
                        (-7 + 8 * K + 3 * K**2) * P041 * Pf1
                        - (27 + 45 * K + 15 * K**2 + K**3) * P042 * Pf1
                        + (
                            (7 - 8 * K - 3 * K**2) * P031
                            + (27 + 45 * K + 15 * K**2 + K**3) * P032
                        )
                        * Pf2
                    )
                    + P011
                    * (
                        (7 - 8 * K - 3 * K**2) * P043 * Pf1
                        - (27 + 45 * K + 15 * K**2 + K**3) * P044 * Pf1
                        + (
                            (-7 + 8 * K + 3 * K**2) * P033
                            + (27 + 45 * K + 15 * K**2 + K**3) * P034
                        )
                        * Pf2
                    )
                )
            )
            * S1**2
        )
        - 2
        * B1**3
        * (
            C1**2
            * (
                -2
                * d1
                * kj**2
                * (6 + 20 * K + 3 * K**2)
                * (
                    P014 * P041 * Pf1
                    + P013 * P042 * Pf1
                    - P012 * P043 * Pf1
                    - P011 * P044 * Pf1
                    - P014 * P031 * Pf2
                    - P013 * P032 * Pf2
                    + P012 * P033 * Pf2
                    + P011 * P034 * Pf2
                )
                + d1**2
                * kj**2
                * (
                    (-(20 + 28 * K + 3 * K**2)) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * kj
                    * (4 + 16 * K + 3 * K**2)
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    - 4
                    * kj**2
                    * K
                    * (4 + 3 * K)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + P012
                    * (
                        (20 + 28 * K + 3 * K**2) * P041 * Pf1
                        - (20 + 28 * K + 3 * K**2) * P031 * Pf2
                        - 2 * kj * (4 + 16 * K + 3 * K**2) * (P044 * Pf1 - P034 * Pf2)
                    )
                )
                - 2
                * K
                * (3 + K)
                * (
                    -2 * P011 * P042 * Pf1
                    + 2 * P011 * P032 * Pf2
                    + 4
                    * kj**2
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + P012
                    * (
                        2 * P041 * Pf1
                        - 2 * P031 * Pf2
                        + kj
                        * (-3 * P043 * Pf1 - P044 * Pf1 + 3 * P033 * Pf2 + P034 * Pf2)
                    )
                    + kj
                    * (
                        P013
                        * ((-P041) * Pf1 + 3 * P042 * Pf1 + (P031 - 3 * P032) * Pf2)
                        + P014
                        * (-3 * P041 * Pf1 + P042 * Pf1 + 3 * P031 * Pf2 - P032 * Pf2)
                        + P011 * (P043 * Pf1 + 3 * P044 * Pf1 - (P033 + 3 * P034) * Pf2)
                    )
                )
            )
            + C1
            * (
                -9 * (4 + 7 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                + 12
                * kj**2
                * K**2
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + P012
                * (
                    9 * (4 + 7 * K + K**2) * P041 * Pf1
                    - 9 * (4 + 7 * K + K**2) * P031 * Pf2
                    - 2
                    * kj
                    * (
                        (12 + 22 * K + 5 * K**2) * P043 * Pf1
                        + 3 * K * P044 * Pf1
                        - 12 * P033 * Pf2
                        - 22 * K * P033 * Pf2
                        - 5 * K**2 * P033 * Pf2
                        - 3 * K * P034 * Pf2
                    )
                )
                + 2
                * kj
                * (
                    -12
                    * (
                        P014 * P041 * Pf1
                        - P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    - 5
                    * K**2
                    * (
                        P014 * P041 * Pf1
                        - P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + K
                    * (
                        P013
                        * (
                            -3 * P041 * Pf1
                            + 22 * P042 * Pf1
                            + 3 * P031 * Pf2
                            - 22 * P032 * Pf2
                        )
                        + P014
                        * (
                            -22 * P041 * Pf1
                            + 3 * P042 * Pf1
                            + 22 * P031 * Pf2
                            - 3 * P032 * Pf2
                        )
                        + P011
                        * (
                            3 * P043 * Pf1
                            + 22 * P044 * Pf1
                            - 3 * P033 * Pf2
                            - 22 * P034 * Pf2
                        )
                    )
                )
            )
            * S1
            + (
                (68 + 51 * K + 4 * K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                + 2
                * d1**2
                * kj**3
                * (4 + 16 * K + 3 * K**2)
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                + 4
                * d1**2
                * kj**4
                * K
                * (4 + 3 * K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + d1
                * kj**2
                * (
                    2 * (6 + 20 * K + 3 * K**2) * P014 * (P041 * Pf1 - P031 * Pf2)
                    + d1 * (20 + 28 * K + 3 * K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 2
                    * (6 + 20 * K + 3 * K**2)
                    * (
                        P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                )
                - 2
                * kj
                * (
                    4
                    * (
                        P014 * P042 * Pf1
                        + P011 * P043 * Pf1
                        - P014 * P032 * Pf2
                        - P011 * P033 * Pf2
                    )
                    - 9
                    * K
                    * (
                        P014 * P041 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + P013
                    * (
                        (-4 + K**2) * P041 * Pf1
                        + 4 * P031 * Pf2
                        + 9 * K * (P042 * Pf1 - P032 * Pf2)
                        + K**2 * (3 * P042 * Pf1 - (P031 + 3 * P032) * Pf2)
                    )
                    + K**2
                    * (
                        P014
                        * (-3 * P041 * Pf1 - P042 * Pf1 + 3 * P031 * Pf2 + P032 * Pf2)
                        + P011
                        * ((-P043) * Pf1 + 3 * P044 * Pf1 + P033 * Pf2 - 3 * P034 * Pf2)
                    )
                )
                + P012
                * (
                    (
                        -(
                            68
                            + 51 * K
                            + 4 * K**2
                            + d1**2 * kj**2 * (20 + 28 * K + 3 * K**2)
                        )
                    )
                    * P041
                    * Pf1
                    + 18 * kj * K * P043 * Pf1
                    + 6 * kj * K**2 * P043 * Pf1
                    + 8 * kj * P044 * Pf1
                    - 2 * kj * K**2 * P044 * Pf1
                    + 68 * P031 * Pf2
                    + 51 * K * P031 * Pf2
                    + 4 * K**2 * P031 * Pf2
                    - 18 * kj * K * P033 * Pf2
                    - 6 * kj * K**2 * P033 * Pf2
                    - 8 * kj * P034 * Pf2
                    + 2 * kj * K**2 * P034 * Pf2
                    - 2
                    * d1
                    * kj**2
                    * (6 + 20 * K + 3 * K**2)
                    * (P043 * Pf1 - P033 * Pf2)
                    + d1**2
                    * kj**2
                    * (
                        (20 + 28 * K + 3 * K**2) * P031 * Pf2
                        + 2 * kj * (4 + 16 * K + 3 * K**2) * (P044 * Pf1 - P034 * Pf2)
                    )
                )
            )
            * S1**2
        )
        + 2
        * B1
        * B2**2
        * (
            2
            * C1**2
            * (
                (-d1)
                * kj**2
                * (12 + 16 * K + 5 * K**2)
                * (
                    P014 * P041 * Pf1
                    + P013 * P042 * Pf1
                    - P011 * P044 * Pf1
                    - P014 * P031 * Pf2
                    - P013 * P032 * Pf2
                    + P011 * P034 * Pf2
                )
                - 3
                * (3 + K)
                * (
                    (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        P014 * P041 * Pf1
                        - P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                )
                + 2
                * d1**2
                * kj**2
                * (2 + K)
                * (
                    (-(2 + K)) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * kj
                    * (1 + K)
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    - 4
                    * kj**2
                    * K
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + P012
                * (
                    (2 + K) * (2 * d1**2 * kj**2 * (2 + K) + 3 * (3 + K)) * P041 * Pf1
                    + d1 * kj**2 * (12 + 16 * K + 5 * K**2) * (P043 * Pf1 - P033 * Pf2)
                    - 3
                    * (3 + K)
                    * ((2 + K) * P031 * Pf2 + kj * K * (P043 * Pf1 - P033 * Pf2))
                    - 2
                    * d1**2
                    * kj**2
                    * (2 + K)
                    * (
                        (2 + K) * P031 * Pf2
                        + 2 * kj * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                    )
                )
            )
            + C1
            * (
                24
                * kj
                * (
                    P014 * P041 * Pf1
                    - P013 * P042 * Pf1
                    + P012 * P043 * Pf1
                    - P011 * P044 * Pf1
                    - P014 * P031 * Pf2
                    + P013 * P032 * Pf2
                    - P012 * P033 * Pf2
                    + P011 * P034 * Pf2
                )
                - 4
                * K
                * (
                    -9 * P011 * P042 * Pf1
                    + 9 * P011 * P032 * Pf2
                    + 6
                    * kj**2
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + kj
                    * (
                        P013
                        * (
                            -3 * P041 * Pf1
                            + 8 * P042 * Pf1
                            + 3 * P031 * Pf2
                            - 8 * P032 * Pf2
                        )
                        + P014
                        * (
                            -8 * P041 * Pf1
                            + 3 * P042 * Pf1
                            + 8 * P031 * Pf2
                            - 3 * P032 * Pf2
                        )
                        + P011
                        * (
                            3 * P043 * Pf1
                            + 8 * P044 * Pf1
                            - 3 * P033 * Pf2
                            - 8 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        9 * P041 * Pf1
                        - 9 * P031 * Pf2
                        + kj
                        * (
                            -8 * P043 * Pf1
                            - 3 * P044 * Pf1
                            + 8 * P033 * Pf2
                            + 3 * P034 * Pf2
                        )
                    )
                )
                + K**2
                * (
                    7 * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 20
                    * kj**2
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + 2
                    * kj
                    * (
                        P014
                        * (5 * P041 * Pf1 - P042 * Pf1 - 5 * P031 * Pf2 + P032 * Pf2)
                        + P013
                        * (P041 * Pf1 - 5 * P042 * Pf1 - P031 * Pf2 + 5 * P032 * Pf2)
                        + P011
                        * ((-P043) * Pf1 - 5 * P044 * Pf1 + P033 * Pf2 + 5 * P034 * Pf2)
                    )
                    + P012
                    * (
                        -7 * P041 * Pf1
                        + 7 * P031 * Pf2
                        + 2
                        * kj
                        * (5 * P043 * Pf1 + P044 * Pf1 - (5 * P033 + P034) * Pf2)
                    )
                )
            )
            * S1
            + 2
            * (
                (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                + 4
                * d1**2
                * kj**3
                * (2 + 3 * K + K**2)
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                + 8
                * d1**2
                * kj**4
                * K
                * (2 + K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + kj**2
                * (2 + K)
                * (
                    2 * d1**2 * (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + d1
                    * (6 + 5 * K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + 8
                    * K
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + P012
                * (
                    (-(2 + K)) * (1 + 2 * d1**2 * kj**2 * (2 + K)) * P041 * Pf1
                    - 9 * kj * K * P043 * Pf1
                    - 3 * kj * K**2 * P043 * Pf1
                    - 4 * kj * P044 * Pf1
                    - 4 * kj * K * P044 * Pf1
                    - kj * K**2 * P044 * Pf1
                    + 2 * P031 * Pf2
                    + K * P031 * Pf2
                    + 9 * kj * K * P033 * Pf2
                    + 3 * kj * K**2 * P033 * Pf2
                    + 4 * kj * P034 * Pf2
                    + 4 * kj * K * P034 * Pf2
                    + kj * K**2 * P034 * Pf2
                    - d1 * kj**2 * (12 + 16 * K + 5 * K**2) * (P043 * Pf1 - P033 * Pf2)
                    + 2
                    * d1**2
                    * kj**2
                    * (2 + K)
                    * (
                        (2 + K) * P031 * Pf2
                        + 2 * kj * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                    )
                )
                + kj
                * (
                    4
                    * (
                        P014 * P042 * Pf1
                        + P011 * P043 * Pf1
                        - P014 * P032 * Pf2
                        - P011 * P033 * Pf2
                    )
                    + P013
                    * (
                        (-((2 + K) ** 2)) * P041 * Pf1
                        + 4 * P031 * Pf2
                        + K * (9 * P042 * Pf1 + 4 * P031 * Pf2 - 9 * P032 * Pf2)
                        + K**2 * (3 * P042 * Pf1 + P031 * Pf2 - 3 * P032 * Pf2)
                    )
                    + K
                    * (
                        P014
                        * (
                            -9 * P041 * Pf1
                            + 4 * P042 * Pf1
                            + 9 * P031 * Pf2
                            - 4 * P032 * Pf2
                        )
                        + P011
                        * (
                            4 * P043 * Pf1
                            + 9 * P044 * Pf1
                            - 4 * P033 * Pf2
                            - 9 * P034 * Pf2
                        )
                    )
                    + K**2
                    * (
                        P014
                        * (-3 * P041 * Pf1 + P042 * Pf1 + 3 * P031 * Pf2 - P032 * Pf2)
                        + P011 * (P043 * Pf1 + 3 * P044 * Pf1 - (P033 + 3 * P034) * Pf2)
                    )
                )
            )
            * S1**2
        )
    )

    N6 = 4 * (
        B1
        * B2
        * (
            C1**2
            * (
                -3 * (24 + 11 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                - 16
                * d1**2
                * kj**3
                * (1 + K)
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                - 16
                * d1**2
                * kj**4
                * (-1 + 2 * K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                - 2
                * kj**2
                * (
                    2 * d1**2 * (5 + 2 * K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (4 + 3 * K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    - 6
                    * K
                    * (3 + K)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + P012
                * (
                    (4 * d1**2 * kj**2 * (5 + 2 * K) + 3 * (24 + 11 * K + K**2))
                    * P041
                    * Pf1
                    - 3 * (24 + 11 * K + K**2) * P031 * Pf2
                    - 2
                    * d1
                    * kj**2
                    * (
                        -3 * (4 + 3 * K) * P043 * Pf1
                        + (2 * d1 * (5 + 2 * K) * P031 + 3 * (4 + 3 * K) * P033) * Pf2
                    )
                    - 16 * d1**2 * kj**3 * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                    - 6
                    * kj
                    * (3 + K)
                    * (
                        (3 + K) * P043 * Pf1
                        + P044 * Pf1
                        - ((3 + K) * P033 + P034) * Pf2
                    )
                )
                + 6
                * kj
                * (3 + K)
                * (
                    P014
                    * (
                        (-(3 + K)) * P041 * Pf1
                        + P042 * Pf1
                        + ((3 + K) * P031 - P032) * Pf2
                    )
                    + P013
                    * (
                        (-P041) * Pf1
                        + (3 + K) * P042 * Pf1
                        + (P031 - (3 + K) * P032) * Pf2
                    )
                    + P011
                    * (
                        P043 * Pf1
                        + (3 + K) * P044 * Pf1
                        - (P033 + (3 + K) * P034) * Pf2
                    )
                )
            )
            + C1
            * (
                (108 + 57 * K + 7 * K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                - 4
                * kj**2
                * (18 + 21 * K + 5 * K**2)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + 2
                * kj
                * (
                    P013
                    * (
                        (18 + 9 * K + K**2) * P041 * Pf1
                        - 6 * (5 + 7 * K + K**2) * P042 * Pf1
                        - ((18 + 9 * K + K**2) * P031 - 6 * (5 + 7 * K + K**2) * P032)
                        * Pf2
                    )
                    + P014
                    * (
                        6 * (5 + 7 * K + K**2) * P041 * Pf1
                        - (18 + 9 * K + K**2) * P042 * Pf1
                        - (6 * (5 + 7 * K + K**2) * P031 - (18 + 9 * K + K**2) * P032)
                        * Pf2
                    )
                    + P011
                    * (
                        (-(18 + 9 * K + K**2)) * P043 * Pf1
                        - 6 * (5 + 7 * K + K**2) * P044 * Pf1
                        + ((18 + 9 * K + K**2) * P033 + 6 * (5 + 7 * K + K**2) * P034)
                        * Pf2
                    )
                )
                + P012
                * (
                    (-(108 + 57 * K + 7 * K**2)) * P041 * Pf1
                    + (108 + 57 * K + 7 * K**2) * P031 * Pf2
                    + 2
                    * kj
                    * (
                        6 * (5 + 7 * K + K**2) * P043 * Pf1
                        + (18 + 9 * K + K**2) * P044 * Pf1
                        - (6 * (5 + 7 * K + K**2) * P033 + (18 + 9 * K + K**2) * P034)
                        * Pf2
                    )
                )
            )
            * S1
            + 2
            * (
                -2 * (1 + 8 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                + 8
                * d1**2
                * kj**3
                * (1 + K)
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                + 8
                * d1**2
                * kj**4
                * (-1 + 2 * K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + kj**2
                * (
                    2 * d1**2 * (5 + 2 * K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (4 + 3 * K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + 4
                    * (1 + 8 * K + K**2)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + kj
                * (
                    P013
                    * (
                        (-(1 + 8 * K + K**2)) * P041 * Pf1
                        + 3 * (3 + K) ** 2 * P042 * Pf1
                        + ((1 + 8 * K + K**2) * P031 - 3 * (3 + K) ** 2 * P032) * Pf2
                    )
                    + P014
                    * (
                        -3 * (3 + K) ** 2 * P041 * Pf1
                        + (1 + 8 * K + K**2) * P042 * Pf1
                        + (3 * (3 + K) ** 2 * P031 - (1 + 8 * K + K**2) * P032) * Pf2
                    )
                    + P011
                    * (
                        (1 + 8 * K + K**2) * P043 * Pf1
                        + 3 * (3 + K) ** 2 * P044 * Pf1
                        - ((1 + 8 * K + K**2) * P033 + 3 * (3 + K) ** 2 * P034) * Pf2
                    )
                )
                + P012
                * (
                    2 * (1 + 8 * K + K**2 - d1**2 * kj**2 * (5 + 2 * K)) * P041 * Pf1
                    - 2 * (1 + 8 * K + K**2) * P031 * Pf2
                    + d1
                    * kj**2
                    * (
                        -3 * (4 + 3 * K) * P043 * Pf1
                        + (2 * d1 * (5 + 2 * K) * P031 + 3 * (4 + 3 * K) * P033) * Pf2
                    )
                    + 8 * d1**2 * kj**3 * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                    + kj
                    * (
                        -3 * (3 + K) ** 2 * P043 * Pf1
                        - (1 + 8 * K + K**2) * P044 * Pf1
                        + (3 * (3 + K) ** 2 * P033 + (1 + 8 * K + K**2) * P034) * Pf2
                    )
                )
            )
            * S1**2
        )
        + B1**2
        * (
            C1**2
            * (
                -2 * (9 + 12 * K + K**2) * P011 * (P042 * Pf1 - P032 * Pf2)
                + 8
                * d1**2
                * kj**3
                * (4 + 3 * K)
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                + 16
                * d1**2
                * kj**4
                * (1 + 3 * K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                - kj
                * (9 + 12 * K + K**2)
                * (
                    P014 * (3 * P041 * Pf1 - P042 * Pf1 - 3 * P031 * Pf2 + P032 * Pf2)
                    + P013 * (P041 * Pf1 - 3 * P042 * Pf1 - P031 * Pf2 + 3 * P032 * Pf2)
                    + P011
                    * ((-P043) * Pf1 - 3 * P044 * Pf1 + P033 * Pf2 + 3 * P034 * Pf2)
                )
                + 2
                * kj**2
                * (
                    2 * d1**2 * (7 + 3 * K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (8 + 5 * K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + 2
                    * (9 + 12 * K + K**2)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + P012
                * (
                    2
                    * (9 + 12 * K + K**2 - 2 * d1**2 * kj**2 * (7 + 3 * K))
                    * P041
                    * Pf1
                    - 2 * (9 + 12 * K + K**2) * P031 * Pf2
                    + 2
                    * d1
                    * kj**2
                    * (
                        -3 * (8 + 5 * K) * P043 * Pf1
                        + (2 * d1 * (7 + 3 * K) * P031 + 3 * (8 + 5 * K) * P033) * Pf2
                    )
                    + 8 * d1**2 * kj**3 * (4 + 3 * K) * (P044 * Pf1 - P034 * Pf2)
                    - kj
                    * (9 + 12 * K + K**2)
                    * (3 * P043 * Pf1 + P044 * Pf1 - (3 * P033 + P034) * Pf2)
                )
            )
            - 3
            * C1
            * (
                (-(36 + 21 * K + K**2)) * P011 * (P042 * Pf1 - P032 * Pf2)
                + 4
                * kj**2
                * K
                * (5 + K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + P012
                * (
                    (36 + 21 * K + K**2) * P041 * Pf1
                    - (36 + 21 * K + K**2) * P031 * Pf2
                    - 2
                    * kj
                    * (
                        (13 + 11 * K + K**2) * P043 * Pf1
                        + (3 + 2 * K) * P044 * Pf1
                        - ((13 + 11 * K + K**2) * P033 + (3 + 2 * K) * P034) * Pf2
                    )
                )
                - 2
                * kj
                * (
                    P014
                    * (
                        (13 + 11 * K + K**2) * P041 * Pf1
                        - (3 + 2 * K) * P042 * Pf1
                        - ((13 + 11 * K + K**2) * P031 - (3 + 2 * K) * P032) * Pf2
                    )
                    + P013
                    * (
                        (3 + 2 * K) * P041 * Pf1
                        - (13 + 11 * K + K**2) * P042 * Pf1
                        + ((-(3 + 2 * K)) * P031 + (13 + 11 * K + K**2) * P032) * Pf2
                    )
                    + P011
                    * (
                        (-(3 + 2 * K)) * P043 * Pf1
                        - (13 + 11 * K + K**2) * P044 * Pf1
                        + ((3 + 2 * K) * P033 + (13 + 11 * K + K**2) * P034) * Pf2
                    )
                )
            )
            * S1
            + (
                (-(106 + 35 * K + K**2)) * P011 * (P042 * Pf1 - P032 * Pf2)
                - 8
                * d1**2
                * kj**3
                * (4 + 3 * K)
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                - 16
                * d1**2
                * kj**4
                * (1 + 3 * K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                - 2
                * kj**2
                * (
                    2 * d1**2 * (7 + 3 * K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (8 + 5 * K)
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    - 2
                    * (-1 + K + 2 * K**2)
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + P012
                * (
                    (106 + 35 * K + K**2 + 4 * d1**2 * kj**2 * (7 + 3 * K)) * P041 * Pf1
                    - (106 + 35 * K + K**2) * P031 * Pf2
                    - 2
                    * d1
                    * kj**2
                    * (
                        -3 * (8 + 5 * K) * P043 * Pf1
                        + (2 * d1 * (7 + 3 * K) * P031 + 3 * (8 + 5 * K) * P033) * Pf2
                    )
                    - 8 * d1**2 * kj**3 * (4 + 3 * K) * (P044 * Pf1 - P034 * Pf2)
                    - kj
                    * (
                        3 * (9 + 12 * K + K**2) * P043 * Pf1
                        - (-17 + 2 * K + K**2) * P044 * Pf1
                        - (3 * (9 + 12 * K + K**2) * P033 - (-17 + 2 * K + K**2) * P034)
                        * Pf2
                    )
                )
                + kj
                * (
                    P014
                    * (
                        -3 * (9 + 12 * K + K**2) * P041 * Pf1
                        - (-17 + 2 * K + K**2) * P042 * Pf1
                        + (3 * (9 + 12 * K + K**2) * P031 + (-17 + 2 * K + K**2) * P032)
                        * Pf2
                    )
                    + P013
                    * (
                        (-17 + 2 * K + K**2) * P041 * Pf1
                        + 3 * (9 + 12 * K + K**2) * P042 * Pf1
                        - ((-17 + 2 * K + K**2) * P031 + 3 * (9 + 12 * K + K**2) * P032)
                        * Pf2
                    )
                    + P011
                    * (
                        (-(-17 + 2 * K + K**2)) * P043 * Pf1
                        + 3 * (9 + 12 * K + K**2) * P044 * Pf1
                        + ((-17 + 2 * K + K**2) * P033 - 3 * (9 + 12 * K + K**2) * P034)
                        * Pf2
                    )
                )
            )
            * S1**2
        )
        + B2**2
        * (
            C1**2
            * (
                -12
                * d1
                * kj**2
                * (2 + K)
                * (
                    P014 * P041 * Pf1
                    + P013 * P042 * Pf1
                    - P011 * P044 * Pf1
                    - P014 * P031 * Pf2
                    - P013 * P032 * Pf2
                    + P011 * P034 * Pf2
                )
                - 9
                * (
                    (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        P014 * P041 * Pf1
                        - P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        + P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                )
                + 4
                * d1**2
                * kj**2
                * (2 + K)
                * (
                    (-P011) * P042 * Pf1
                    + P011 * P032 * Pf2
                    - 2
                    * kj
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    - 4
                    * kj**2
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + P012
                * (
                    (9 + 4 * d1**2 * kj**2) * (2 + K) * P041 * Pf1
                    + 12 * d1 * kj**2 * (2 + K) * (P043 * Pf1 - P033 * Pf2)
                    - 4
                    * d1**2
                    * kj**2
                    * (2 + K)
                    * (2 * kj * P044 * Pf1 + P031 * Pf2 - 2 * kj * P034 * Pf2)
                    - 9 * ((2 + K) * P031 * Pf2 + kj * K * (P043 * Pf1 - P033 * Pf2))
                )
            )
            - 6
            * C1
            * (
                -4
                * kj
                * (
                    P014 * P041 * Pf1
                    - P013 * P042 * Pf1
                    + P012 * P043 * Pf1
                    - P011 * P044 * Pf1
                    - P014 * P031 * Pf2
                    + P013 * P032 * Pf2
                    - P012 * P033 * Pf2
                    + P011 * P034 * Pf2
                )
                + K
                * (
                    -2 * P011 * P042 * Pf1
                    + 2 * P011 * P032 * Pf2
                    + 4
                    * kj**2
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + P012
                    * (
                        2 * P041 * Pf1
                        - 2 * P031 * Pf2
                        + kj
                        * (-2 * P043 * Pf1 - P044 * Pf1 + 2 * P033 * Pf2 + P034 * Pf2)
                    )
                    + kj
                    * (
                        P013
                        * ((-P041) * Pf1 + 2 * P042 * Pf1 + (P031 - 2 * P032) * Pf2)
                        + P014
                        * (-2 * P041 * Pf1 + P042 * Pf1 + 2 * P031 * Pf2 - P032 * Pf2)
                        + P011 * (P043 * Pf1 + 2 * P044 * Pf1 - (P033 + 2 * P034) * Pf2)
                    )
                )
            )
            * S1
            + (
                (2 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                + 8
                * d1**2
                * kj**3
                * (2 + K)
                * (
                    P013 * P041 * Pf1
                    - P014 * P042 * Pf1
                    - P011 * P043 * Pf1
                    - P013 * P031 * Pf2
                    + P014 * P032 * Pf2
                    + P011 * P033 * Pf2
                )
                + 16
                * d1**2
                * kj**4
                * (2 + K)
                * (
                    P014 * P043 * Pf1
                    - P013 * P044 * Pf1
                    - P014 * P033 * Pf2
                    + P013 * P034 * Pf2
                )
                + P012
                * (
                    (-(1 + 4 * d1**2 * kj**2)) * (2 + K) * P041 * Pf1
                    - 9 * kj * K * P043 * Pf1
                    - 8 * kj * P044 * Pf1
                    - 4 * kj * K * P044 * Pf1
                    + 2 * P031 * Pf2
                    + K * P031 * Pf2
                    + 9 * kj * K * P033 * Pf2
                    + 8 * kj * P034 * Pf2
                    + 4 * kj * K * P034 * Pf2
                    - 12 * d1 * kj**2 * (2 + K) * (P043 * Pf1 - P033 * Pf2)
                    + 4
                    * d1**2
                    * kj**2
                    * (2 + K)
                    * (2 * kj * P044 * Pf1 + P031 * Pf2 - 2 * kj * P034 * Pf2)
                )
                + 4
                * kj**2
                * (2 + K)
                * (
                    d1**2 * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + 4
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                )
                + kj
                * (
                    P013
                    * (
                        -4 * (2 + K) * P041 * Pf1
                        + 9 * K * P042 * Pf1
                        + 8 * P031 * Pf2
                        + 4 * K * P031 * Pf2
                        - 9 * K * P032 * Pf2
                    )
                    + 8
                    * (
                        P014 * P042 * Pf1
                        + P011 * P043 * Pf1
                        - P014 * P032 * Pf2
                        - P011 * P033 * Pf2
                    )
                    + K
                    * (
                        P014
                        * (
                            -9 * P041 * Pf1
                            + 4 * P042 * Pf1
                            + 9 * P031 * Pf2
                            - 4 * P032 * Pf2
                        )
                        + P011
                        * (
                            4 * P043 * Pf1
                            + 9 * P044 * Pf1
                            - 4 * P033 * Pf2
                            - 9 * P034 * Pf2
                        )
                    )
                )
            )
            * S1**2
        )
    )

    N7 = (
        -3
        * (C1 - S1)
        * (
            2
            * B1
            * (
                2
                * C1
                * (
                    2 * (3 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * d1**2
                    * kj**3
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    - 8
                    * d1**2
                    * kj**4
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + kj
                    * (3 + K)
                    * (
                        P014
                        * (3 * P041 * Pf1 - P042 * Pf1 - 3 * P031 * Pf2 + P032 * Pf2)
                        + P013
                        * (P041 * Pf1 - 3 * P042 * Pf1 - P031 * Pf2 + 3 * P032 * Pf2)
                        + P011
                        * ((-P043) * Pf1 - 3 * P044 * Pf1 + P033 * Pf2 + 3 * P034 * Pf2)
                    )
                    - 2
                    * kj**2
                    * (
                        d1**2 * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 3
                        * d1
                        * (
                            P014 * P041 * Pf1
                            + P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            - P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                        + 2
                        * (3 + K)
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        2 * (-3 + d1**2 * kj**2 - K) * P041 * Pf1
                        + 2 * (3 + K) * P031 * Pf2
                        - 2
                        * d1
                        * kj**2
                        * (-3 * P043 * Pf1 + d1 * P031 * Pf2 + 3 * P033 * Pf2)
                        + 4 * d1**2 * kj**3 * ((-P044) * Pf1 + P034 * Pf2)
                        + kj
                        * (3 + K)
                        * (3 * P043 * Pf1 + P044 * Pf1 - (3 * P033 + P034) * Pf2)
                    )
                )
                + (
                    -3 * (8 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    - 8
                    * d1**2
                    * kj**3
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    - 16
                    * d1**2
                    * kj**4
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    - 4
                    * kj**2
                    * (
                        d1**2 * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 3
                        * d1
                        * (
                            P014 * P041 * Pf1
                            + P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            - P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                        + 3
                        * K
                        * (
                            (-P014) * P043 * Pf1
                            + P013 * P044 * Pf1
                            + P014 * P033 * Pf2
                            - P013 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        (4 * d1**2 * kj**2 + 3 * (8 + K)) * P041 * Pf1
                        - 3 * (8 + K) * P031 * Pf2
                        - 4
                        * d1
                        * kj**2
                        * (-3 * P043 * Pf1 + d1 * P031 * Pf2 + 3 * P033 * Pf2)
                        - 8 * d1**2 * kj**3 * (P044 * Pf1 - P034 * Pf2)
                        - 6
                        * kj
                        * (
                            (3 + K) * P043 * Pf1
                            + P044 * Pf1
                            - ((3 + K) * P033 + P034) * Pf2
                        )
                    )
                    - 6
                    * kj
                    * (
                        P014
                        * (
                            (3 + K) * P041 * Pf1
                            - P042 * Pf1
                            - 3 * P031 * Pf2
                            - K * P031 * Pf2
                            + P032 * Pf2
                        )
                        + P013
                        * (
                            P041 * Pf1
                            - (3 + K) * P042 * Pf1
                            + (-P031 + (3 + K) * P032) * Pf2
                        )
                        + P011
                        * (
                            (-P043) * Pf1
                            - (3 + K) * P044 * Pf1
                            + (P033 + (3 + K) * P034) * Pf2
                        )
                    )
                )
                * S1
            )
            + B2
            * (
                C1
                * (
                    3 * (8 + K) * P011 * (P042 * Pf1 - P032 * Pf2)
                    + 16
                    * d1**2
                    * kj**3
                    * (
                        P013 * P041 * Pf1
                        - P014 * P042 * Pf1
                        - P011 * P043 * Pf1
                        - P013 * P031 * Pf2
                        + P014 * P032 * Pf2
                        + P011 * P033 * Pf2
                    )
                    + 32
                    * d1**2
                    * kj**4
                    * (
                        P014 * P043 * Pf1
                        - P013 * P044 * Pf1
                        - P014 * P033 * Pf2
                        + P013 * P034 * Pf2
                    )
                    + 4
                    * kj**2
                    * (
                        2 * d1**2 * P011 * (P042 * Pf1 - P032 * Pf2)
                        + 6
                        * d1
                        * (
                            P014 * P041 * Pf1
                            + P013 * P042 * Pf1
                            - P011 * P044 * Pf1
                            - P014 * P031 * Pf2
                            - P013 * P032 * Pf2
                            + P011 * P034 * Pf2
                        )
                        + 3
                        * K
                        * (
                            (-P014) * P043 * Pf1
                            + P013 * P044 * Pf1
                            + P014 * P033 * Pf2
                            - P013 * P034 * Pf2
                        )
                    )
                    + P012
                    * (
                        (-(8 * d1**2 * kj**2 + 3 * (8 + K))) * P041 * Pf1
                        + 3 * (8 + K) * P031 * Pf2
                        + 8
                        * d1
                        * kj**2
                        * (-3 * P043 * Pf1 + d1 * P031 * Pf2 + 3 * P033 * Pf2)
                        + 16 * d1**2 * kj**3 * (P044 * Pf1 - P034 * Pf2)
                        + 6
                        * kj
                        * (
                            (3 + K) * P043 * Pf1
                            + P044 * Pf1
                            - ((3 + K) * P033 + P034) * Pf2
                        )
                    )
                    + 6
                    * kj
                    * (
                        P014
                        * (
                            (3 + K) * P041 * Pf1
                            - P042 * Pf1
                            - 3 * P031 * Pf2
                            - K * P031 * Pf2
                            + P032 * Pf2
                        )
                        + P013
                        * (
                            P041 * Pf1
                            - (3 + K) * P042 * Pf1
                            + (-P031 + (3 + K) * P032) * Pf2
                        )
                        + P011
                        * (
                            (-P043) * Pf1
                            - (3 + K) * P044 * Pf1
                            + (P033 + (3 + K) * P034) * Pf2
                        )
                    )
                )
                + (
                    24
                    * d1
                    * kj**2
                    * (
                        P014 * P041 * Pf1
                        + P013 * P042 * Pf1
                        - P012 * P043 * Pf1
                        - P011 * P044 * Pf1
                        - P014 * P031 * Pf2
                        - P013 * P032 * Pf2
                        + P012 * P033 * Pf2
                        + P011 * P034 * Pf2
                    )
                    + 6
                    * kj
                    * (
                        P011 * P043 * Pf1
                        - 3 * P012 * P043 * Pf1
                        + 3 * P011 * P044 * Pf1
                        - P012 * P044 * Pf1
                        - P011 * P033 * Pf2
                        + 3 * P012 * P033 * Pf2
                        - 3 * P011 * P034 * Pf2
                        + P012 * P034 * Pf2
                        + P014
                        * (
                            -3 * P041 * Pf1
                            + P042 * Pf1
                            + 8 * kj * P043 * Pf1
                            + 3 * P031 * Pf2
                            - P032 * Pf2
                            - 8 * kj * P033 * Pf2
                        )
                        + P013
                        * (
                            (-P041) * Pf1
                            + 3 * P042 * Pf1
                            - 8 * kj * P044 * Pf1
                            + P031 * Pf2
                            - 3 * P032 * Pf2
                            + 8 * kj * P034 * Pf2
                        )
                    )
                    - 8
                    * d1**2
                    * kj**2
                    * (
                        (-P011) * P042 * Pf1
                        + P011 * P032 * Pf2
                        - 2
                        * kj
                        * (
                            P013 * P041 * Pf1
                            - P014 * P042 * Pf1
                            - P011 * P043 * Pf1
                            - P013 * P031 * Pf2
                            + P014 * P032 * Pf2
                            + P011 * P033 * Pf2
                        )
                        + P012
                        * (
                            P041 * Pf1
                            - 2 * kj * P044 * Pf1
                            - P031 * Pf2
                            + 2 * kj * P034 * Pf2
                        )
                        - 4
                        * kj**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                    )
                    + K
                    * (
                        -5 * P011 * P042 * Pf1
                        + 5 * P011 * P032 * Pf2
                        + 4
                        * kj**2
                        * (
                            P014 * P043 * Pf1
                            - P013 * P044 * Pf1
                            - P014 * P033 * Pf2
                            + P013 * P034 * Pf2
                        )
                        + kj
                        * (
                            P013
                            * (
                                -4 * P041 * Pf1
                                + 6 * P042 * Pf1
                                + 4 * P031 * Pf2
                                - 6 * P032 * Pf2
                            )
                            + P014
                            * (
                                -6 * P041 * Pf1
                                + 4 * P042 * Pf1
                                + 6 * P031 * Pf2
                                - 4 * P032 * Pf2
                            )
                            + 2
                            * P011
                            * (
                                2 * P043 * Pf1
                                + 3 * P044 * Pf1
                                - 2 * P033 * Pf2
                                - 3 * P034 * Pf2
                            )
                        )
                        + P012
                        * (
                            5 * P041 * Pf1
                            - 5 * P031 * Pf2
                            + kj
                            * (
                                -6 * P043 * Pf1
                                - 4 * P044 * Pf1
                                + 6 * P033 * Pf2
                                + 4 * P034 * Pf2
                            )
                        )
                    )
                )
                * S1
            )
        )
    )

    N8 = (
        9
        * (
            -2 * P011 * P042 * Pf1
            + 2 * P011 * P032 * Pf2
            + 4
            * kj**2
            * (
                P014 * P043 * Pf1
                - P013 * P044 * Pf1
                - P014 * P033 * Pf2
                + P013 * P034 * Pf2
            )
            + P012
            * (
                2 * P041 * Pf1
                - 2 * P031 * Pf2
                + kj * (-3 * P043 * Pf1 - P044 * Pf1 + 3 * P033 * Pf2 + P034 * Pf2)
            )
            + kj
            * (
                P013 * ((-P041) * Pf1 + 3 * P042 * Pf1 + (P031 - 3 * P032) * Pf2)
                + P014 * (-3 * P041 * Pf1 + P042 * Pf1 + 3 * P031 * Pf2 - P032 * Pf2)
                + P011 * (P043 * Pf1 + 3 * P044 * Pf1 - (P033 + 3 * P034) * Pf2)
            )
        )
        * (C1 - S1) ** 2
    )

    return [N1, N2, N3, N4, N5, N6, N7, N8]


def get_num2_2visc(P0, Pzsf, d1, kj, B1, B2, K):

    P011, P012, P013, P014 = P0[0, :]
    P021, P022, P023, P024 = P0[1, :]
    P031, P032, P033, P034 = P0[2, :]
    P041, P042, P043, P044 = P0[3, :]

    # Calculating hyperbolic functions
    C1 = cosh(kj * d1)
    S1 = sinh(kj * d1)
    Pf1 = Pzsf[0]
    Pf2 = Pzsf[1]

    N1 = (
        64
        * B1**4
        * B2
        * K**3
        * (
            P022 * P041 * Pf1
            - P021 * P042 * Pf1
            - P022 * P031 * Pf2
            + P021 * P032 * Pf2
        )
        * (
            -2 * B1 * B2 * C1 * S1
            + B2**2 * (C1**2 * (1 + d1**2 * kj**2) - d1**2 * kj**2 * S1**2)
            + B1**2 * ((-(C1**2)) * d1**2 * kj**2 + (1 + d1**2 * kj**2) * S1**2)
        )
    )

    N2 = (
        -64
        * B1**3
        * K**2
        * (
            B1**3
            * (
                P022 * P041 * Pf1
                - P021 * P042 * Pf1
                - P022 * P031 * Pf2
                + P021 * P032 * Pf2
            )
            * (C1**2 * d1**2 * kj**2 - (1 + d1**2 * kj**2) * S1**2)
            + B1
            * B2**2
            * (
                C1**2
                * (
                    (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + d1**2 * kj**2 * (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        P024 * P041 * Pf1
                        - P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + P022
                    * (
                        (-(1 + d1**2 * kj**2)) * (2 + K) * P041 * Pf1
                        + (2 + K) * P031 * Pf2
                        + d1**2 * kj**2 * (2 + K) * P031 * Pf2
                        + kj * K * (P043 * Pf1 - P033 * Pf2)
                    )
                )
                + C1
                * (7 + 3 * K)
                * (
                    P022 * P041 * Pf1
                    - P021 * P042 * Pf1
                    - P022 * P031 * Pf2
                    + P021 * P032 * Pf2
                )
                * S1
                + kj
                * (
                    d1**2
                    * kj
                    * (2 + K)
                    * (
                        P022 * P041 * Pf1
                        - P021 * P042 * Pf1
                        - P022 * P031 * Pf2
                        + P021 * P032 * Pf2
                    )
                    + K
                    * (
                        P024 * P041 * Pf1
                        - P023 * P042 * Pf1
                        + P022 * P043 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        - P022 * P033 * Pf2
                        + P021 * P034 * Pf2
                    )
                )
                * S1**2
            )
            + B1**2
            * B2
            * (
                C1**2
                * d1
                * kj**2
                * (
                    K
                    * (
                        (-P024) * P041 * Pf1
                        - P023 * P042 * Pf1
                        + P022 * P043 * Pf1
                        + P021 * P044 * Pf1
                        + P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        - P022 * P033 * Pf2
                        - P021 * P034 * Pf2
                    )
                    + d1
                    * (
                        (-(3 + 2 * K)) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            (-P023) * P041 * Pf1
                            + P024 * P042 * Pf1
                            + P021 * P043 * Pf1
                            + P023 * P031 * Pf2
                            - P024 * P032 * Pf2
                            - P021 * P033 * Pf2
                        )
                        + P022
                        * (
                            (3 + 2 * K) * P041 * Pf1
                            - (3 + 2 * K) * P031 * Pf2
                            + kj * K * ((-P044) * Pf1 + P034 * Pf2)
                        )
                    )
                )
                + C1
                * (
                    (-(3 + K)) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        (-P024) * P041 * Pf1
                        + P023 * P042 * Pf1
                        + P021 * P044 * Pf1
                        + P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        - P021 * P034 * Pf2
                    )
                    + P022
                    * (
                        (3 + K) * P041 * Pf1
                        - (3 + K) * P031 * Pf2
                        + kj * K * ((-P043) * Pf1 + P033 * Pf2)
                    )
                )
                * S1
                + (
                    2 * (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + d1
                    * kj**2
                    * K
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (3 + 2 * K) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                    )
                    + P022
                    * (
                        (-(2 * (2 + K) + d1**2 * kj**2 * (3 + 2 * K))) * P041 * Pf1
                        + 2 * (2 + K) * P031 * Pf2
                        + d1 * kj**2 * K * ((-P043) * Pf1 + P033 * Pf2)
                        + d1**2
                        * kj**2
                        * (
                            (3 + 2 * K) * P031 * Pf2
                            + kj * K * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B2**3
            * (
                C1**2
                * (
                    (3 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + d1
                    * kj**2
                    * K
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                    )
                    + P022
                    * (
                        (-(3 + K + d1**2 * kj**2 * (2 + K))) * P041 * Pf1
                        + (3 + K) * P031 * Pf2
                        + d1 * kj**2 * K * ((-P043) * Pf1 + P033 * Pf2)
                        + d1**2
                        * kj**2
                        * ((2 + K) * P031 * Pf2 + kj * K * (P044 * Pf1 - P034 * Pf2))
                    )
                )
                + C1
                * kj
                * K
                * (
                    (-P024) * P041 * Pf1
                    + P023 * P042 * Pf1
                    - P022 * P043 * Pf1
                    + P021 * P044 * Pf1
                    + P024 * P031 * Pf2
                    - P023 * P032 * Pf2
                    + P022 * P033 * Pf2
                    - P021 * P034 * Pf2
                )
                * S1
                + d1
                * kj**2
                * (
                    K
                    * (
                        (-P024) * P041 * Pf1
                        - P023 * P042 * Pf1
                        + P022 * P043 * Pf1
                        + P021 * P044 * Pf1
                        + P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        - P022 * P033 * Pf2
                        - P021 * P034 * Pf2
                    )
                    + d1
                    * (
                        (-(2 + K)) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            (-P023) * P041 * Pf1
                            + P024 * P042 * Pf1
                            + P021 * P043 * Pf1
                            + P023 * P031 * Pf2
                            - P024 * P032 * Pf2
                            - P021 * P033 * Pf2
                        )
                        + P022
                        * (
                            (2 + K) * P041 * Pf1
                            - (2 + K) * P031 * Pf2
                            + kj * K * ((-P044) * Pf1 + P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
        )
    )

    N3 = (
        16
        * B1**2
        * K
        * (
            B1**2
            * B2
            * (
                C1**2
                * (
                    2
                    * d1
                    * kj**2
                    * K
                    * (7 + 3 * K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P022 * P043 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P022 * P033 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (12 + 14 * K + 5 * K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 6
                        * kj
                        * K
                        * (2 + K)
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        + 4
                        * kj**2
                        * K**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                        + P022
                        * (
                            (-(12 + 14 * K + 5 * K**2)) * P041 * Pf1
                            + (12 + 14 * K + 5 * K**2) * P031 * Pf2
                            + 6 * kj * K * (2 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                    + K
                    * (
                        (-(8 + K)) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 4
                        * kj**2
                        * K
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                        + P022
                        * (
                            (8 + K) * P041 * Pf1
                            - (8 + K) * P031 * Pf2
                            - 2
                            * kj
                            * (
                                (3 + K) * P043 * Pf1
                                + P044 * Pf1
                                - ((3 + K) * P033 + P034) * Pf2
                            )
                        )
                        - 2
                        * kj
                        * (
                            P024
                            * (
                                (3 + K) * P041 * Pf1
                                - P042 * Pf1
                                - 3 * P031 * Pf2
                                - K * P031 * Pf2
                                + P032 * Pf2
                            )
                            + P023
                            * (
                                P041 * Pf1
                                - (3 + K) * P042 * Pf1
                                + (-P031 + (3 + K) * P032) * Pf2
                            )
                            + P021
                            * (
                                (-P043) * Pf1
                                - (3 + K) * P044 * Pf1
                                + (P033 + (3 + K) * P034) * Pf2
                            )
                        )
                    )
                )
                - 2
                * C1
                * (7 + 3 * K)
                * (
                    (-(3 + K)) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        (-P024) * P041 * Pf1
                        + P023 * P042 * Pf1
                        + P021 * P044 * Pf1
                        + P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        - P021 * P034 * Pf2
                    )
                    + P022
                    * (
                        (3 + K) * P041 * Pf1
                        - (3 + K) * P031 * Pf2
                        + kj * K * ((-P043) * Pf1 + P033 * Pf2)
                    )
                )
                * S1
                + (
                    -2
                    * d1
                    * kj**2
                    * K
                    * (7 + 3 * K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    - 2
                    * (
                        (9 + 16 * K + 3 * K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (3 + K)
                        * (
                            P024 * P041 * Pf1
                            - P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            + P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                    )
                    + d1**2
                    * kj**2
                    * (
                        (-(12 + 14 * K + 5 * K**2)) * P021 * (P042 * Pf1 - P032 * Pf2)
                        - 6
                        * kj
                        * K
                        * (2 + K)
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        - 4
                        * kj**2
                        * K**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        (
                            18
                            + 32 * K
                            + 6 * K**2
                            + d1**2 * kj**2 * (12 + 14 * K + 5 * K**2)
                        )
                        * P041
                        * Pf1
                        - 2 * (9 + 16 * K + 3 * K**2) * P031 * Pf2
                        - 2 * kj * K * (3 + K) * (P043 * Pf1 - P033 * Pf2)
                        - 6 * d1**2 * kj**3 * K * (2 + K) * (P044 * Pf1 - P034 * Pf2)
                        - d1
                        * kj**2
                        * (
                            12 * d1 * P031 * Pf2
                            + K**2
                            * (-6 * P043 * Pf1 + 5 * d1 * P031 * Pf2 + 6 * P033 * Pf2)
                            - 14 * K * (P043 * Pf1 - (d1 * P031 + P033) * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + 2
            * B1
            * B2**2
            * (
                2
                * C1**2
                * (
                    (-d1)
                    * kj**2
                    * K
                    * (2 + K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    - d1**2
                    * kj**2
                    * (2 + K)
                    * (
                        (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                    )
                    - (3 + K)
                    * (
                        (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P024 * P041 * Pf1
                            - P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            + P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        (2 + K) * (3 + K + d1**2 * kj**2 * (2 + K)) * P041 * Pf1
                        - (6 + 5 * K + K**2) * P031 * Pf2
                        - kj * K * (3 + K) * (P043 * Pf1 - P033 * Pf2)
                        - d1**2 * kj**3 * K * (2 + K) * (P044 * Pf1 - P034 * Pf2)
                        - d1
                        * kj**2
                        * (2 + K)
                        * (
                            2 * d1 * P031 * Pf2
                            + K * ((-P043) * Pf1 + (d1 * P031 + P033) * Pf2)
                        )
                    )
                )
                + C1
                * (
                    3 * (4 + 7 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * kj**2
                    * K**2
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + P022
                    * (
                        -3 * (4 + 7 * K + K**2) * P041 * Pf1
                        + 3 * (4 + 7 * K + K**2) * P031 * Pf2
                        + 2
                        * kj
                        * K
                        * (
                            (2 + K) * P043 * Pf1
                            + P044 * Pf1
                            - ((2 + K) * P033 + P034) * Pf2
                        )
                    )
                    - 2
                    * kj
                    * K
                    * (
                        P024
                        * (
                            (-(2 + K)) * P041 * Pf1
                            + P042 * Pf1
                            + ((2 + K) * P031 - P032) * Pf2
                        )
                        + P023
                        * (
                            (-P041) * Pf1
                            + (2 + K) * P042 * Pf1
                            + (P031 - (2 + K) * P032) * Pf2
                        )
                        + P021
                        * (
                            P043 * Pf1
                            + (2 + K) * P044 * Pf1
                            - (P033 + (2 + K) * P034) * Pf2
                        )
                    )
                )
                * S1
                + 2
                * kj
                * (
                    (-K)
                    * (3 + K)
                    * (
                        P024 * P041 * Pf1
                        - P023 * P042 * Pf1
                        + P022 * P043 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        - P022 * P033 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1
                    * kj
                    * K
                    * (2 + K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P022 * P043 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P022 * P033 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1**2
                    * kj
                    * (2 + K)
                    * (
                        (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        + P022
                        * (
                            (-(2 + K)) * P041 * Pf1
                            + (2 + K) * P031 * Pf2
                            + kj * K * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + 2
            * B1**3
            * (
                C1**2
                * d1
                * kj**2
                * (
                    2
                    * K
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P022 * P043 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P022 * P033 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1
                    * (
                        (4 + 5 * K) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 2
                        * kj
                        * K
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        + P022
                        * (
                            (-(4 + 5 * K)) * P041 * Pf1
                            + (4 + 5 * K) * P031 * Pf2
                            + 2 * kj * K * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                + C1
                * K
                * (
                    3 * P021 * (P042 * Pf1 - P032 * Pf2)
                    + P022
                    * (
                        -3 * P041 * Pf1
                        + 2 * kj * P043 * Pf1
                        + 3 * P031 * Pf2
                        - 2 * kj * P033 * Pf2
                    )
                    + 2
                    * kj
                    * (
                        P024 * P041 * Pf1
                        - P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                )
                * S1
                + (
                    -4 * (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * d1
                    * kj**2
                    * K
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (-(4 + 5 * K)) * P021 * (P042 * Pf1 - P032 * Pf2)
                        - 2
                        * kj
                        * K
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                    )
                    + P022
                    * (
                        (4 * (2 + K) + d1**2 * kj**2 * (4 + 5 * K)) * P041 * Pf1
                        - 4 * (2 + K) * P031 * Pf2
                        + 2 * d1 * kj**2 * K * (P043 * Pf1 - P033 * Pf2)
                        - d1**2
                        * kj**2
                        * (
                            (4 + 5 * K) * P031 * Pf2
                            + 2 * kj * K * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B2**3
            * (
                C1**2
                * (
                    (-(9 + 12 * K + K**2)) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * d1
                    * kj**2
                    * K
                    * (5 + K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (-(4 + 8 * K + K**2)) * P021 * (P042 * Pf1 - P032 * Pf2)
                        - 2
                        * kj
                        * K
                        * (4 + K)
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        - 4
                        * kj**2
                        * K**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        (9 + 12 * K + K**2 + d1**2 * kj**2 * (4 + 8 * K + K**2))
                        * P041
                        * Pf1
                        - (9 + 12 * K + K**2) * P031 * Pf2
                        + 2 * d1 * kj**2 * K * (5 + K) * (P043 * Pf1 - P033 * Pf2)
                        - d1**2
                        * kj**2
                        * (
                            (4 + 8 * K + K**2) * P031 * Pf2
                            + 2 * kj * K * (4 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                + 2
                * C1
                * kj
                * K
                * (5 + K)
                * (
                    P024 * P041 * Pf1
                    - P023 * P042 * Pf1
                    + P022 * P043 * Pf1
                    - P021 * P044 * Pf1
                    - P024 * P031 * Pf2
                    + P023 * P032 * Pf2
                    - P022 * P033 * Pf2
                    + P021 * P034 * Pf2
                )
                * S1
                + (
                    P021 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * kj
                    * K
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    + 2
                    * d1**2
                    * kj**3
                    * K
                    * (4 + K)
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    + 4
                    * d1**2
                    * kj**4
                    * K**2
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + kj**2
                    * (
                        d1**2 * (4 + 8 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 2
                        * d1
                        * K
                        * (5 + K)
                        * (
                            P024 * P041 * Pf1
                            + P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            - P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                        + 4
                        * K**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        (-(1 + d1**2 * kj**2 * (4 + 8 * K + K**2))) * P041 * Pf1
                        - 2 * kj * K * P044 * Pf1
                        + P031 * Pf2
                        + 2 * kj * K * P034 * Pf2
                        - 2 * d1 * kj**2 * K * (5 + K) * (P043 * Pf1 - P033 * Pf2)
                        + d1**2
                        * kj**2
                        * (
                            (4 + 8 * K + K**2) * P031 * Pf2
                            + 2 * kj * K * (4 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
        )
    )

    N4 = (
        16
        * B1
        * (
            B1
            * B2**2
            * (
                C1**2
                * (
                    -2
                    * d1
                    * kj**2
                    * K
                    * (10 + 7 * K + K**2)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    - (9 + 12 * K + K**2)
                    * (
                        (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + kj
                        * K
                        * (
                            P024 * P041 * Pf1
                            - P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            + P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                    )
                    + d1**2
                    * kj**2
                    * (2 + K)
                    * (
                        (-(4 + 8 * K + K**2)) * P021 * (P042 * Pf1 - P032 * Pf2)
                        - 2
                        * kj
                        * K
                        * (4 + K)
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        - 4
                        * kj**2
                        * K**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        (2 + K)
                        * (9 + 12 * K + K**2 + d1**2 * kj**2 * (4 + 8 * K + K**2))
                        * P041
                        * Pf1
                        - (18 + 33 * K + 14 * K**2 + K**3) * P031 * Pf2
                        - kj * K * (9 + 12 * K + K**2) * (P043 * Pf1 - P033 * Pf2)
                        - 2
                        * d1**2
                        * kj**3
                        * K
                        * (8 + 6 * K + K**2)
                        * (P044 * Pf1 - P034 * Pf2)
                        - d1
                        * kj**2
                        * (2 + K)
                        * (
                            4 * d1 * P031 * Pf2
                            + K**2
                            * (-2 * P043 * Pf1 + d1 * P031 * Pf2 + 2 * P033 * Pf2)
                            + 2
                            * K
                            * (-5 * P043 * Pf1 + 4 * d1 * P031 * Pf2 + 5 * P033 * Pf2)
                        )
                    )
                )
                + C1
                * K
                * (
                    (36 + 21 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * kj**2
                    * K
                    * (5 + K)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + P022
                    * (
                        (-(36 + 21 * K + K**2)) * P041 * Pf1
                        + (36 + 21 * K + K**2) * P031 * Pf2
                        + 2
                        * kj
                        * (
                            (10 + 7 * K + K**2) * P043 * Pf1
                            + (3 + 2 * K) * P044 * Pf1
                            - ((10 + 7 * K + K**2) * P033 + (3 + 2 * K) * P034) * Pf2
                        )
                    )
                    + 2
                    * kj
                    * (
                        P024
                        * (
                            (10 + 7 * K + K**2) * P041 * Pf1
                            - (3 + 2 * K) * P042 * Pf1
                            - ((10 + 7 * K + K**2) * P031 - (3 + 2 * K) * P032) * Pf2
                        )
                        + P023
                        * (
                            (3 + 2 * K) * P041 * Pf1
                            - (10 + 7 * K + K**2) * P042 * Pf1
                            + ((-(3 + 2 * K)) * P031 + (10 + 7 * K + K**2) * P032) * Pf2
                        )
                        + P021
                        * (
                            (-(3 + 2 * K)) * P043 * Pf1
                            - (10 + 7 * K + K**2) * P044 * Pf1
                            + ((3 + 2 * K) * P033 + (10 + 7 * K + K**2) * P034) * Pf2
                        )
                    )
                )
                * S1
                + (
                    (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 2
                    * d1**2
                    * kj**3
                    * K
                    * (8 + 6 * K + K**2)
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    + 4
                    * d1**2
                    * kj**4
                    * K**2
                    * (2 + K)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + kj**2
                    * (2 + K)
                    * (
                        d1**2 * (4 + 8 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 2
                        * d1
                        * K
                        * (5 + K)
                        * (
                            P024 * P041 * Pf1
                            + P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            - P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                        + 4
                        * K**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                    )
                    + kj
                    * K
                    * (
                        P024
                        * (
                            (-(9 + 12 * K + K**2)) * P041 * Pf1
                            + 2 * (2 + K) * P042 * Pf1
                            + ((9 + 12 * K + K**2) * P031 - 2 * (2 + K) * P032) * Pf2
                        )
                        + P023
                        * (
                            -2 * (2 + K) * P041 * Pf1
                            + (9 + 12 * K + K**2) * P042 * Pf1
                            + (2 * (2 + K) * P031 - (9 + 12 * K + K**2) * P032) * Pf2
                        )
                        + P021
                        * (
                            2 * (2 + K) * P043 * Pf1
                            + (9 + 12 * K + K**2) * P044 * Pf1
                            - (2 * (2 + K) * P033 + (9 + 12 * K + K**2) * P034) * Pf2
                        )
                    )
                    + P022
                    * (
                        (-(2 + K))
                        * (1 + d1**2 * kj**2 * (4 + 8 * K + K**2))
                        * P041
                        * Pf1
                        + (2 + K) * P031 * Pf2
                        + 2
                        * d1**2
                        * kj**3
                        * K
                        * (8 + 6 * K + K**2)
                        * (P044 * Pf1 - P034 * Pf2)
                        + kj
                        * K
                        * (
                            (-(9 + 12 * K + K**2)) * P043 * Pf1
                            - 2 * (2 + K) * P044 * Pf1
                            + ((9 + 12 * K + K**2) * P033 + 2 * (2 + K) * P034) * Pf2
                        )
                        + d1
                        * kj**2
                        * (2 + K)
                        * (
                            4 * d1 * P031 * Pf2
                            + K**2
                            * (-2 * P043 * Pf1 + d1 * P031 * Pf2 + 2 * P033 * Pf2)
                            + 2
                            * K
                            * (-5 * P043 * Pf1 + 4 * d1 * P031 * Pf2 + 5 * P033 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B1**2
            * B2
            * (
                C1**2
                * (
                    d1
                    * kj**2
                    * K
                    * (16 + 11 * K + 2 * K**2)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P022 * P043 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P022 * P033 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (4 + 4 * K + 5 * K**2 + K**3) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 2
                        * kj
                        * K
                        * (6 + 4 * K + K**2)
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        + 4
                        * kj**2
                        * K**2
                        * (3 + K)
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                        + P022
                        * (
                            (-(4 + 4 * K + 5 * K**2 + K**3)) * P041 * Pf1
                            + (4 + 4 * K + 5 * K**2 + K**3) * P031 * Pf2
                            + 2
                            * kj
                            * K
                            * (6 + 4 * K + K**2)
                            * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                    + K
                    * (3 + K)
                    * (
                        (-(8 + K)) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 4
                        * kj**2
                        * K
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                        + P022
                        * (
                            (8 + K) * P041 * Pf1
                            - (8 + K) * P031 * Pf2
                            - 2
                            * kj
                            * (
                                (3 + K) * P043 * Pf1
                                + P044 * Pf1
                                - ((3 + K) * P033 + P034) * Pf2
                            )
                        )
                        - 2
                        * kj
                        * (
                            P024
                            * (
                                (3 + K) * P041 * Pf1
                                - P042 * Pf1
                                - 3 * P031 * Pf2
                                - K * P031 * Pf2
                                + P032 * Pf2
                            )
                            + P023
                            * (
                                P041 * Pf1
                                - (3 + K) * P042 * Pf1
                                + (-P031 + (3 + K) * P032) * Pf2
                            )
                            + P021
                            * (
                                (-P043) * Pf1
                                - (3 + K) * P044 * Pf1
                                + (P033 + (3 + K) * P034) * Pf2
                            )
                        )
                    )
                )
                + C1
                * (
                    3
                    * (12 + 25 * K + 10 * K**2 + K**3)
                    * P021
                    * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * kj**2
                    * K**2
                    * (3 + K)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + P022
                    * (
                        -3 * (12 + 25 * K + 10 * K**2 + K**3) * P041 * Pf1
                        + 3 * (12 + 25 * K + 10 * K**2 + K**3) * P031 * Pf2
                        + kj
                        * K
                        * (
                            (14 + 29 * K + 4 * K**2) * P043 * Pf1
                            + 2 * (3 + K) * P044 * Pf1
                            - ((14 + 29 * K + 4 * K**2) * P033 + 2 * (3 + K) * P034)
                            * Pf2
                        )
                    )
                    + kj
                    * K
                    * (
                        P024
                        * (
                            (14 + 29 * K + 4 * K**2) * P041 * Pf1
                            - 2 * (3 + K) * P042 * Pf1
                            - ((14 + 29 * K + 4 * K**2) * P031 - 2 * (3 + K) * P032)
                            * Pf2
                        )
                        + P023
                        * (
                            2 * (3 + K) * P041 * Pf1
                            - (14 + 29 * K + 4 * K**2) * P042 * Pf1
                            + (-2 * (3 + K) * P031 + (14 + 29 * K + 4 * K**2) * P032)
                            * Pf2
                        )
                        + P021
                        * (
                            -2 * (3 + K) * P043 * Pf1
                            - (14 + 29 * K + 4 * K**2) * P044 * Pf1
                            + (2 * (3 + K) * P033 + (14 + 29 * K + 4 * K**2) * P034)
                            * Pf2
                        )
                    )
                )
                * S1
                + (
                    (-(4 + 35 * K + 24 * K**2 + 2 * K**3))
                    * P021
                    * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * d1**2
                    * kj**3
                    * K
                    * (6 + 4 * K + K**2)
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    - 4
                    * d1**2
                    * kj**4
                    * K**2
                    * (3 + K)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + d1
                    * kj**2
                    * (
                        4 * d1 * P021 * ((-P042) * Pf1 + P032 * Pf2)
                        + K**2
                        * (
                            -11 * P024 * P041 * Pf1
                            - 5 * d1 * P021 * P042 * Pf1
                            - 11 * P023 * P042 * Pf1
                            + 11 * P021 * P044 * Pf1
                            + 11 * P024 * P031 * Pf2
                            + 5 * d1 * P021 * P032 * Pf2
                            + 11 * P023 * P032 * Pf2
                            - 11 * P021 * P034 * Pf2
                        )
                        + K**3
                        * (
                            -2 * P024 * P041 * Pf1
                            - d1 * P021 * P042 * Pf1
                            - 2 * P023 * P042 * Pf1
                            + 2 * P021 * P044 * Pf1
                            + 2 * P024 * P031 * Pf2
                            + d1 * P021 * P032 * Pf2
                            + 2 * P023 * P032 * Pf2
                            - 2 * P021 * P034 * Pf2
                        )
                        - 4
                        * K
                        * (
                            4 * P024 * P041 * Pf1
                            + d1 * P021 * P042 * Pf1
                            + 4 * P023 * P042 * Pf1
                            - 4 * P021 * P044 * Pf1
                            - 4 * P024 * P031 * Pf2
                            - d1 * P021 * P032 * Pf2
                            - 4 * P023 * P032 * Pf2
                            + 4 * P021 * P034 * Pf2
                        )
                    )
                    + kj
                    * K
                    * (
                        P024
                        * (
                            -2 * (3 + K) ** 2 * P041 * Pf1
                            + (-2 + K) * P042 * Pf1
                            + (2 * (3 + K) ** 2 * P031 - (-2 + K) * P032) * Pf2
                        )
                        + P023
                        * (
                            (-(-2 + K)) * P041 * Pf1
                            + 2 * (3 + K) ** 2 * P042 * Pf1
                            + ((-2 + K) * P031 - 2 * (3 + K) ** 2 * P032) * Pf2
                        )
                        + P021
                        * (
                            (-2 + K) * P043 * Pf1
                            + 2 * (3 + K) ** 2 * P044 * Pf1
                            - ((-2 + K) * P033 + 2 * (3 + K) ** 2 * P034) * Pf2
                        )
                    )
                    + P022
                    * (
                        (
                            4
                            + 35 * K
                            + 24 * K**2
                            + 2 * K**3
                            + d1**2 * kj**2 * (4 + 4 * K + 5 * K**2 + K**3)
                        )
                        * P041
                        * Pf1
                        - (4 + 35 * K + 24 * K**2 + 2 * K**3) * P031 * Pf2
                        - 2
                        * d1**2
                        * kj**3
                        * K
                        * (6 + 4 * K + K**2)
                        * (P044 * Pf1 - P034 * Pf2)
                        - kj
                        * K
                        * (
                            2 * (3 + K) ** 2 * P043 * Pf1
                            + (-2 + K) * P044 * Pf1
                            - (2 * (3 + K) ** 2 * P033 + (-2 + K) * P034) * Pf2
                        )
                        - d1
                        * kj**2
                        * (
                            4 * d1 * P031 * Pf2
                            + K**3
                            * (-2 * P043 * Pf1 + d1 * P031 * Pf2 + 2 * P033 * Pf2)
                            + 4
                            * K
                            * (-4 * P043 * Pf1 + d1 * P031 * Pf2 + 4 * P033 * Pf2)
                            + K**2
                            * (-11 * P043 * Pf1 + 5 * d1 * P031 * Pf2 + 11 * P033 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B2**3
            * K
            * (
                C1**2
                * (
                    -3 * (3 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - d1
                    * kj**2
                    * (6 + 5 * K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + 2
                    * d1**2
                    * kj**2
                    * (
                        (-(2 + K)) * P021 * (P042 * Pf1 - P032 * Pf2)
                        - 2
                        * kj
                        * (1 + K)
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        - 4
                        * kj**2
                        * K
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        (2 * d1**2 * kj**2 * (2 + K) + 3 * (3 + K)) * P041 * Pf1
                        - 3 * (3 + K) * P031 * Pf2
                        + d1 * kj**2 * (6 + 5 * K) * (P043 * Pf1 - P033 * Pf2)
                        - 2
                        * d1**2
                        * kj**2
                        * (
                            (2 + K) * P031 * Pf2
                            + 2 * kj * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                + C1
                * kj
                * (6 + 5 * K)
                * (
                    P024 * P041 * Pf1
                    - P023 * P042 * Pf1
                    + P022 * P043 * Pf1
                    - P021 * P044 * Pf1
                    - P024 * P031 * Pf2
                    + P023 * P032 * Pf2
                    - P022 * P033 * Pf2
                    + P021 * P034 * Pf2
                )
                * S1
                + (
                    P021 * (P042 * Pf1 - P032 * Pf2)
                    + 4
                    * d1**2
                    * kj**3
                    * (1 + K)
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    - kj
                    * (2 + K)
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    + 8
                    * d1**2
                    * kj**4
                    * K
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + kj**2
                    * (
                        2 * d1**2 * (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + d1
                        * (6 + 5 * K)
                        * (
                            P024 * P041 * Pf1
                            + P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            - P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                        + 8
                        * K
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        (-(1 + 2 * d1**2 * kj**2 * (2 + K))) * P041 * Pf1
                        - 2 * kj * P044 * Pf1
                        - kj * K * P044 * Pf1
                        + P031 * Pf2
                        + 2 * kj * P034 * Pf2
                        + kj * K * P034 * Pf2
                        - d1 * kj**2 * (6 + 5 * K) * (P043 * Pf1 - P033 * Pf2)
                        + 2
                        * d1**2
                        * kj**2
                        * (
                            (2 + K) * P031 * Pf2
                            + 2 * kj * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
            + B1**3
            * (
                C1**2
                * (
                    2
                    * d1
                    * kj**2
                    * K
                    * (5 + 4 * K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P022 * P043 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P022 * P033 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + d1**2
                    * kj**2
                    * (
                        (4 + 20 * K + 7 * K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 8
                        * kj
                        * K
                        * (1 + K)
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        + 4
                        * kj**2
                        * K**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                        + P022
                        * (
                            (-(4 + 20 * K + 7 * K**2)) * P041 * Pf1
                            + (4 + 20 * K + 7 * K**2) * P031 * Pf2
                            + 8 * kj * K * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                    + K**2
                    * (
                        -2 * P021 * P042 * Pf1
                        + 2 * P021 * P032 * Pf2
                        + 4
                        * kj**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                        + P022
                        * (
                            2 * P041 * Pf1
                            - 2 * P031 * Pf2
                            + kj
                            * (
                                -3 * P043 * Pf1
                                - P044 * Pf1
                                + 3 * P033 * Pf2
                                + P034 * Pf2
                            )
                        )
                        + kj
                        * (
                            P023
                            * ((-P041) * Pf1 + 3 * P042 * Pf1 + (P031 - 3 * P032) * Pf2)
                            + P024
                            * (
                                -3 * P041 * Pf1
                                + P042 * Pf1
                                + 3 * P031 * Pf2
                                - P032 * Pf2
                            )
                            + P021
                            * (P043 * Pf1 + 3 * P044 * Pf1 - (P033 + 3 * P034) * Pf2)
                        )
                    )
                )
                + C1
                * K
                * (7 + 3 * K)
                * (
                    3 * P021 * (P042 * Pf1 - P032 * Pf2)
                    + P022
                    * (
                        -3 * P041 * Pf1
                        + 2 * kj * P043 * Pf1
                        + 3 * P031 * Pf2
                        - 2 * kj * P033 * Pf2
                    )
                    + 2
                    * kj
                    * (
                        P024 * P041 * Pf1
                        - P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                )
                * S1
                + (
                    (-(16 + 33 * K + 6 * K**2)) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 8
                    * d1**2
                    * kj**3
                    * K
                    * (1 + K)
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    - 4
                    * d1**2
                    * kj**4
                    * K**2
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + d1
                    * kj**2
                    * (
                        4 * d1 * P021 * ((-P042) * Pf1 + P032 * Pf2)
                        + K**2
                        * (
                            -8 * P024 * P041 * Pf1
                            - 7 * d1 * P021 * P042 * Pf1
                            - 8 * P023 * P042 * Pf1
                            + 8 * P021 * P044 * Pf1
                            + 8 * P024 * P031 * Pf2
                            + 7 * d1 * P021 * P032 * Pf2
                            + 8 * P023 * P032 * Pf2
                            - 8 * P021 * P034 * Pf2
                        )
                        - 10
                        * K
                        * (
                            P024 * P041 * Pf1
                            + 2 * d1 * P021 * P042 * Pf1
                            + P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            - 2 * d1 * P021 * P032 * Pf2
                            - P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                    )
                    + kj
                    * K
                    * (
                        2
                        * (
                            P024 * P042 * Pf1
                            + P021 * P043 * Pf1
                            - P024 * P032 * Pf2
                            - P021 * P033 * Pf2
                        )
                        + P023
                        * (
                            (-2 + K) * P041 * Pf1
                            + 2 * P031 * Pf2
                            + K * (3 * P042 * Pf1 - (P031 + 3 * P032) * Pf2)
                        )
                        + K
                        * (
                            P024
                            * (
                                -3 * P041 * Pf1
                                - P042 * Pf1
                                + 3 * P031 * Pf2
                                + P032 * Pf2
                            )
                            + P021
                            * (
                                (-P043) * Pf1
                                + 3 * P044 * Pf1
                                + P033 * Pf2
                                - 3 * P034 * Pf2
                            )
                        )
                    )
                    + P022
                    * (
                        (
                            16
                            + 33 * K
                            + 6 * K**2
                            + d1**2 * kj**2 * (4 + 20 * K + 7 * K**2)
                        )
                        * P041
                        * Pf1
                        - 3 * kj * K**2 * P043 * Pf1
                        - 2 * kj * K * P044 * Pf1
                        + kj * K**2 * P044 * Pf1
                        - 16 * P031 * Pf2
                        - 33 * K * P031 * Pf2
                        - 6 * K**2 * P031 * Pf2
                        + 3 * kj * K**2 * P033 * Pf2
                        + 2 * kj * K * P034 * Pf2
                        - kj * K**2 * P034 * Pf2
                        + 2 * d1 * kj**2 * K * (5 + 4 * K) * (P043 * Pf1 - P033 * Pf2)
                        - d1**2
                        * kj**2
                        * (
                            (4 + 20 * K + 7 * K**2) * P031 * Pf2
                            + 8 * kj * K * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                        )
                    )
                )
                * S1**2
            )
        )
    )

    N5 = 4 * (
        B2**3
        * K
        * (
            C1**2
            * (
                -9 * P021 * P042 * Pf1
                + 9 * P021 * P032 * Pf2
                - 12
                * d1
                * kj**2
                * (
                    P024 * P041 * Pf1
                    + P023 * P042 * Pf1
                    - P021 * P044 * Pf1
                    - P024 * P031 * Pf2
                    - P023 * P032 * Pf2
                    + P021 * P034 * Pf2
                )
                + P022
                * (
                    (9 + 4 * d1**2 * kj**2) * P041 * Pf1
                    - 9 * P031 * Pf2
                    + 12 * d1 * kj**2 * (P043 * Pf1 - P033 * Pf2)
                    - 4
                    * d1**2
                    * kj**2
                    * (2 * kj * P044 * Pf1 + P031 * Pf2 - 2 * kj * P034 * Pf2)
                )
                + 4
                * d1**2
                * kj**2
                * (
                    (-P021) * P042 * Pf1
                    + P021 * P032 * Pf2
                    - 2
                    * kj
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    - 4
                    * kj**2
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
            )
            + 12
            * C1
            * kj
            * (
                P024 * P041 * Pf1
                - P023 * P042 * Pf1
                + P022 * P043 * Pf1
                - P021 * P044 * Pf1
                - P024 * P031 * Pf2
                + P023 * P032 * Pf2
                - P022 * P033 * Pf2
                + P021 * P034 * Pf2
            )
            * S1
            + (
                P021 * (P042 * Pf1 - P032 * Pf2)
                - 4
                * kj
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                + 8
                * d1**2
                * kj**3
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                + 16
                * d1**2
                * kj**4
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + P022
                * (
                    -4 * kj * P044 * Pf1
                    - P041 * (Pf1 + 4 * d1**2 * kj**2 * Pf1)
                    + P031 * Pf2
                    + 4 * kj * P034 * Pf2
                    - 12 * d1 * kj**2 * (P043 * Pf1 - P033 * Pf2)
                    + 4
                    * d1**2
                    * kj**2
                    * (2 * kj * P044 * Pf1 + P031 * Pf2 - 2 * kj * P034 * Pf2)
                )
                + 4
                * kj**2
                * (
                    d1**2 * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + 4
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
            )
            * S1**2
        )
        + B1**2
        * B2
        * (
            C1**2
            * (
                (-(72 + 105 * K + 20 * K**2 + K**3)) * P021 * (P042 * Pf1 - P032 * Pf2)
                + 4
                * d1**2
                * kj**3
                * (4 - 2 * K + K**2)
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                + 8
                * d1**2
                * kj**4
                * K
                * (6 + K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + 2
                * kj**2
                * (
                    d1**2 * (-4 - 10 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 2
                    * d1
                    * (6 - K + 2 * K**2)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + 2
                    * K
                    * (9 + 12 * K + K**2)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + P022
                * (
                    (
                        72
                        + 105 * K
                        + 20 * K**2
                        + K**3
                        - 2 * d1**2 * kj**2 * (-4 - 10 * K + K**2)
                    )
                    * P041
                    * Pf1
                    - (72 + 105 * K + 20 * K**2 + K**3) * P031 * Pf2
                    + 2
                    * d1
                    * kj**2
                    * (
                        -2 * (6 - K + 2 * K**2) * P043 * Pf1
                        + (
                            d1 * (-4 - 10 * K + K**2) * P031
                            + 2 * (6 - K + 2 * K**2) * P033
                        )
                        * Pf2
                    )
                    + 4 * d1**2 * kj**3 * (4 - 2 * K + K**2) * (P044 * Pf1 - P034 * Pf2)
                    - 2
                    * kj
                    * (9 + 12 * K + K**2)
                    * (
                        (3 + K) * P043 * Pf1
                        + P044 * Pf1
                        - ((3 + K) * P033 + P034) * Pf2
                    )
                )
                + 2
                * kj
                * (9 + 12 * K + K**2)
                * (
                    P024
                    * (
                        (-(3 + K)) * P041 * Pf1
                        + P042 * Pf1
                        + ((3 + K) * P031 - P032) * Pf2
                    )
                    + P023
                    * (
                        (-P041) * Pf1
                        + (3 + K) * P042 * Pf1
                        + (P031 - (3 + K) * P032) * Pf2
                    )
                    + P021
                    * (
                        P043 * Pf1
                        + (3 + K) * P044 * Pf1
                        - (P033 + (3 + K) * P034) * Pf2
                    )
                )
            )
            - 2
            * C1
            * (
                (-(108 + 99 * K + 24 * K**2 + K**3)) * P021 * (P042 * Pf1 - P032 * Pf2)
                + 4
                * kj**2
                * K
                * (15 + 8 * K + K**2)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + P022
                * (
                    (108 + 99 * K + 24 * K**2 + K**3) * P041 * Pf1
                    - (108 + 99 * K + 24 * K**2 + K**3) * P031 * Pf2
                    - 2
                    * kj
                    * (
                        (3 + 40 * K + 17 * K**2 + K**3) * P043 * Pf1
                        + (9 + 9 * K + 2 * K**2) * P044 * Pf1
                        - (
                            (3 + 40 * K + 17 * K**2 + K**3) * P033
                            + (9 + 9 * K + 2 * K**2) * P034
                        )
                        * Pf2
                    )
                )
                - 2
                * kj
                * (
                    P024
                    * (
                        (3 + 40 * K + 17 * K**2 + K**3) * P041 * Pf1
                        - (9 + 9 * K + 2 * K**2) * P042 * Pf1
                        - (
                            (3 + 40 * K + 17 * K**2 + K**3) * P031
                            - (9 + 9 * K + 2 * K**2) * P032
                        )
                        * Pf2
                    )
                    + P023
                    * (
                        (9 + 9 * K + 2 * K**2) * P041 * Pf1
                        - (3 + 40 * K + 17 * K**2 + K**3) * P042 * Pf1
                        - (
                            (9 + 9 * K + 2 * K**2) * P031
                            - (3 + 40 * K + 17 * K**2 + K**3) * P032
                        )
                        * Pf2
                    )
                    + P021
                    * (
                        (-(9 + 9 * K + 2 * K**2)) * P043 * Pf1
                        - (3 + 40 * K + 17 * K**2 + K**3) * P044 * Pf1
                        + (
                            (9 + 9 * K + 2 * K**2) * P033
                            + (3 + 40 * K + 17 * K**2 + K**3) * P034
                        )
                        * Pf2
                    )
                )
            )
            * S1
            + (
                (-(16 + 101 * K + 32 * K**2 + K**3)) * P021 * (P042 * Pf1 - P032 * Pf2)
                - 4
                * d1**2
                * kj**3
                * (4 - 2 * K + K**2)
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                - 8
                * d1**2
                * kj**4
                * K
                * (6 + K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + kj**2
                * (
                    -2 * d1**2 * (-4 - 10 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * d1
                    * (6 - K + 2 * K**2)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + 4
                    * K
                    * (1 + 8 * K + K**2)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + P022
                * (
                    (
                        16
                        + 101 * K
                        + 32 * K**2
                        + K**3
                        + 2 * d1**2 * kj**2 * (-4 - 10 * K + K**2)
                    )
                    * P041
                    * Pf1
                    - (16 + 101 * K + 32 * K**2 + K**3) * P031 * Pf2
                    - 2
                    * d1
                    * kj**2
                    * (
                        -2 * (6 - K + 2 * K**2) * P043 * Pf1
                        + (
                            d1 * (-4 - 10 * K + K**2) * P031
                            + 2 * (6 - K + 2 * K**2) * P033
                        )
                        * Pf2
                    )
                    - 4 * d1**2 * kj**3 * (4 - 2 * K + K**2) * (P044 * Pf1 - P034 * Pf2)
                    - 2
                    * kj
                    * (
                        (27 + 45 * K + 15 * K**2 + K**3) * P043 * Pf1
                        + (-7 + 8 * K + 3 * K**2) * P044 * Pf1
                        - (
                            (27 + 45 * K + 15 * K**2 + K**3) * P033
                            + (-7 + 8 * K + 3 * K**2) * P034
                        )
                        * Pf2
                    )
                )
                - 2
                * kj
                * (
                    P024
                    * (
                        (27 + 45 * K + 15 * K**2 + K**3) * P041 * Pf1
                        - (-7 + 8 * K + 3 * K**2) * P042 * Pf1
                        - (
                            (27 + 45 * K + 15 * K**2 + K**3) * P031
                            + (7 - 8 * K - 3 * K**2) * P032
                        )
                        * Pf2
                    )
                    + P023
                    * (
                        (-7 + 8 * K + 3 * K**2) * P041 * Pf1
                        - (27 + 45 * K + 15 * K**2 + K**3) * P042 * Pf1
                        + (
                            (7 - 8 * K - 3 * K**2) * P031
                            + (27 + 45 * K + 15 * K**2 + K**3) * P032
                        )
                        * Pf2
                    )
                    + P021
                    * (
                        (7 - 8 * K - 3 * K**2) * P043 * Pf1
                        - (27 + 45 * K + 15 * K**2 + K**3) * P044 * Pf1
                        + (
                            (-7 + 8 * K + 3 * K**2) * P033
                            + (27 + 45 * K + 15 * K**2 + K**3) * P034
                        )
                        * Pf2
                    )
                )
            )
            * S1**2
        )
        - 2
        * B1**3
        * (
            C1**2
            * (
                -2
                * d1
                * kj**2
                * (6 + 20 * K + 3 * K**2)
                * (
                    P024 * P041 * Pf1
                    + P023 * P042 * Pf1
                    - P022 * P043 * Pf1
                    - P021 * P044 * Pf1
                    - P024 * P031 * Pf2
                    - P023 * P032 * Pf2
                    + P022 * P033 * Pf2
                    + P021 * P034 * Pf2
                )
                + d1**2
                * kj**2
                * (
                    (-(20 + 28 * K + 3 * K**2)) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * kj
                    * (4 + 16 * K + 3 * K**2)
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    - 4
                    * kj**2
                    * K
                    * (4 + 3 * K)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + P022
                    * (
                        (20 + 28 * K + 3 * K**2) * P041 * Pf1
                        - (20 + 28 * K + 3 * K**2) * P031 * Pf2
                        - 2 * kj * (4 + 16 * K + 3 * K**2) * (P044 * Pf1 - P034 * Pf2)
                    )
                )
                - 2
                * K
                * (3 + K)
                * (
                    -2 * P021 * P042 * Pf1
                    + 2 * P021 * P032 * Pf2
                    + 4
                    * kj**2
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + P022
                    * (
                        2 * P041 * Pf1
                        - 2 * P031 * Pf2
                        + kj
                        * (-3 * P043 * Pf1 - P044 * Pf1 + 3 * P033 * Pf2 + P034 * Pf2)
                    )
                    + kj
                    * (
                        P023
                        * ((-P041) * Pf1 + 3 * P042 * Pf1 + (P031 - 3 * P032) * Pf2)
                        + P024
                        * (-3 * P041 * Pf1 + P042 * Pf1 + 3 * P031 * Pf2 - P032 * Pf2)
                        + P021 * (P043 * Pf1 + 3 * P044 * Pf1 - (P033 + 3 * P034) * Pf2)
                    )
                )
            )
            + C1
            * (
                -9 * (4 + 7 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                + 12
                * kj**2
                * K**2
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + P022
                * (
                    9 * (4 + 7 * K + K**2) * P041 * Pf1
                    - 9 * (4 + 7 * K + K**2) * P031 * Pf2
                    - 2
                    * kj
                    * (
                        (12 + 22 * K + 5 * K**2) * P043 * Pf1
                        + 3 * K * P044 * Pf1
                        - 12 * P033 * Pf2
                        - 22 * K * P033 * Pf2
                        - 5 * K**2 * P033 * Pf2
                        - 3 * K * P034 * Pf2
                    )
                )
                + 2
                * kj
                * (
                    -12
                    * (
                        P024 * P041 * Pf1
                        - P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    - 5
                    * K**2
                    * (
                        P024 * P041 * Pf1
                        - P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + K
                    * (
                        P023
                        * (
                            -3 * P041 * Pf1
                            + 22 * P042 * Pf1
                            + 3 * P031 * Pf2
                            - 22 * P032 * Pf2
                        )
                        + P024
                        * (
                            -22 * P041 * Pf1
                            + 3 * P042 * Pf1
                            + 22 * P031 * Pf2
                            - 3 * P032 * Pf2
                        )
                        + P021
                        * (
                            3 * P043 * Pf1
                            + 22 * P044 * Pf1
                            - 3 * P033 * Pf2
                            - 22 * P034 * Pf2
                        )
                    )
                )
            )
            * S1
            + (
                (68 + 51 * K + 4 * K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                + 2
                * d1**2
                * kj**3
                * (4 + 16 * K + 3 * K**2)
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                + 4
                * d1**2
                * kj**4
                * K
                * (4 + 3 * K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + d1
                * kj**2
                * (
                    2 * (6 + 20 * K + 3 * K**2) * P024 * (P041 * Pf1 - P031 * Pf2)
                    + d1 * (20 + 28 * K + 3 * K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 2
                    * (6 + 20 * K + 3 * K**2)
                    * (
                        P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                )
                - 2
                * kj
                * (
                    4
                    * (
                        P024 * P042 * Pf1
                        + P021 * P043 * Pf1
                        - P024 * P032 * Pf2
                        - P021 * P033 * Pf2
                    )
                    - 9
                    * K
                    * (
                        P024 * P041 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + P023
                    * (
                        (-4 + K**2) * P041 * Pf1
                        + 4 * P031 * Pf2
                        + 9 * K * (P042 * Pf1 - P032 * Pf2)
                        + K**2 * (3 * P042 * Pf1 - (P031 + 3 * P032) * Pf2)
                    )
                    + K**2
                    * (
                        P024
                        * (-3 * P041 * Pf1 - P042 * Pf1 + 3 * P031 * Pf2 + P032 * Pf2)
                        + P021
                        * ((-P043) * Pf1 + 3 * P044 * Pf1 + P033 * Pf2 - 3 * P034 * Pf2)
                    )
                )
                + P022
                * (
                    (
                        -(
                            68
                            + 51 * K
                            + 4 * K**2
                            + d1**2 * kj**2 * (20 + 28 * K + 3 * K**2)
                        )
                    )
                    * P041
                    * Pf1
                    + 18 * kj * K * P043 * Pf1
                    + 6 * kj * K**2 * P043 * Pf1
                    + 8 * kj * P044 * Pf1
                    - 2 * kj * K**2 * P044 * Pf1
                    + 68 * P031 * Pf2
                    + 51 * K * P031 * Pf2
                    + 4 * K**2 * P031 * Pf2
                    - 18 * kj * K * P033 * Pf2
                    - 6 * kj * K**2 * P033 * Pf2
                    - 8 * kj * P034 * Pf2
                    + 2 * kj * K**2 * P034 * Pf2
                    - 2
                    * d1
                    * kj**2
                    * (6 + 20 * K + 3 * K**2)
                    * (P043 * Pf1 - P033 * Pf2)
                    + d1**2
                    * kj**2
                    * (
                        (20 + 28 * K + 3 * K**2) * P031 * Pf2
                        + 2 * kj * (4 + 16 * K + 3 * K**2) * (P044 * Pf1 - P034 * Pf2)
                    )
                )
            )
            * S1**2
        )
        + 2
        * B1
        * B2**2
        * (
            2
            * C1**2
            * (
                (-d1)
                * kj**2
                * (12 + 16 * K + 5 * K**2)
                * (
                    P024 * P041 * Pf1
                    + P023 * P042 * Pf1
                    - P021 * P044 * Pf1
                    - P024 * P031 * Pf2
                    - P023 * P032 * Pf2
                    + P021 * P034 * Pf2
                )
                - 3
                * (3 + K)
                * (
                    (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        P024 * P041 * Pf1
                        - P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                )
                + 2
                * d1**2
                * kj**2
                * (2 + K)
                * (
                    (-(2 + K)) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 2
                    * kj
                    * (1 + K)
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    - 4
                    * kj**2
                    * K
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + P022
                * (
                    (2 + K) * (2 * d1**2 * kj**2 * (2 + K) + 3 * (3 + K)) * P041 * Pf1
                    + d1 * kj**2 * (12 + 16 * K + 5 * K**2) * (P043 * Pf1 - P033 * Pf2)
                    - 3
                    * (3 + K)
                    * ((2 + K) * P031 * Pf2 + kj * K * (P043 * Pf1 - P033 * Pf2))
                    - 2
                    * d1**2
                    * kj**2
                    * (2 + K)
                    * (
                        (2 + K) * P031 * Pf2
                        + 2 * kj * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                    )
                )
            )
            + C1
            * (
                24
                * kj
                * (
                    P024 * P041 * Pf1
                    - P023 * P042 * Pf1
                    + P022 * P043 * Pf1
                    - P021 * P044 * Pf1
                    - P024 * P031 * Pf2
                    + P023 * P032 * Pf2
                    - P022 * P033 * Pf2
                    + P021 * P034 * Pf2
                )
                - 4
                * K
                * (
                    -9 * P021 * P042 * Pf1
                    + 9 * P021 * P032 * Pf2
                    + 6
                    * kj**2
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + kj
                    * (
                        P023
                        * (
                            -3 * P041 * Pf1
                            + 8 * P042 * Pf1
                            + 3 * P031 * Pf2
                            - 8 * P032 * Pf2
                        )
                        + P024
                        * (
                            -8 * P041 * Pf1
                            + 3 * P042 * Pf1
                            + 8 * P031 * Pf2
                            - 3 * P032 * Pf2
                        )
                        + P021
                        * (
                            3 * P043 * Pf1
                            + 8 * P044 * Pf1
                            - 3 * P033 * Pf2
                            - 8 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        9 * P041 * Pf1
                        - 9 * P031 * Pf2
                        + kj
                        * (
                            -8 * P043 * Pf1
                            - 3 * P044 * Pf1
                            + 8 * P033 * Pf2
                            + 3 * P034 * Pf2
                        )
                    )
                )
                + K**2
                * (
                    7 * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 20
                    * kj**2
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + 2
                    * kj
                    * (
                        P024
                        * (5 * P041 * Pf1 - P042 * Pf1 - 5 * P031 * Pf2 + P032 * Pf2)
                        + P023
                        * (P041 * Pf1 - 5 * P042 * Pf1 - P031 * Pf2 + 5 * P032 * Pf2)
                        + P021
                        * ((-P043) * Pf1 - 5 * P044 * Pf1 + P033 * Pf2 + 5 * P034 * Pf2)
                    )
                    + P022
                    * (
                        -7 * P041 * Pf1
                        + 7 * P031 * Pf2
                        + 2
                        * kj
                        * (5 * P043 * Pf1 + P044 * Pf1 - (5 * P033 + P034) * Pf2)
                    )
                )
            )
            * S1
            + 2
            * (
                (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                + 4
                * d1**2
                * kj**3
                * (2 + 3 * K + K**2)
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                + 8
                * d1**2
                * kj**4
                * K
                * (2 + K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + kj**2
                * (2 + K)
                * (
                    2 * d1**2 * (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + d1
                    * (6 + 5 * K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + 8
                    * K
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + P022
                * (
                    (-(2 + K)) * (1 + 2 * d1**2 * kj**2 * (2 + K)) * P041 * Pf1
                    - 9 * kj * K * P043 * Pf1
                    - 3 * kj * K**2 * P043 * Pf1
                    - 4 * kj * P044 * Pf1
                    - 4 * kj * K * P044 * Pf1
                    - kj * K**2 * P044 * Pf1
                    + 2 * P031 * Pf2
                    + K * P031 * Pf2
                    + 9 * kj * K * P033 * Pf2
                    + 3 * kj * K**2 * P033 * Pf2
                    + 4 * kj * P034 * Pf2
                    + 4 * kj * K * P034 * Pf2
                    + kj * K**2 * P034 * Pf2
                    - d1 * kj**2 * (12 + 16 * K + 5 * K**2) * (P043 * Pf1 - P033 * Pf2)
                    + 2
                    * d1**2
                    * kj**2
                    * (2 + K)
                    * (
                        (2 + K) * P031 * Pf2
                        + 2 * kj * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                    )
                )
                + kj
                * (
                    4
                    * (
                        P024 * P042 * Pf1
                        + P021 * P043 * Pf1
                        - P024 * P032 * Pf2
                        - P021 * P033 * Pf2
                    )
                    + P023
                    * (
                        (-((2 + K) ** 2)) * P041 * Pf1
                        + 4 * P031 * Pf2
                        + K * (9 * P042 * Pf1 + 4 * P031 * Pf2 - 9 * P032 * Pf2)
                        + K**2 * (3 * P042 * Pf1 + P031 * Pf2 - 3 * P032 * Pf2)
                    )
                    + K
                    * (
                        P024
                        * (
                            -9 * P041 * Pf1
                            + 4 * P042 * Pf1
                            + 9 * P031 * Pf2
                            - 4 * P032 * Pf2
                        )
                        + P021
                        * (
                            4 * P043 * Pf1
                            + 9 * P044 * Pf1
                            - 4 * P033 * Pf2
                            - 9 * P034 * Pf2
                        )
                    )
                    + K**2
                    * (
                        P024
                        * (-3 * P041 * Pf1 + P042 * Pf1 + 3 * P031 * Pf2 - P032 * Pf2)
                        + P021 * (P043 * Pf1 + 3 * P044 * Pf1 - (P033 + 3 * P034) * Pf2)
                    )
                )
            )
            * S1**2
        )
    )

    N6 = 4 * (
        B1
        * B2
        * (
            C1**2
            * (
                -3 * (24 + 11 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                - 16
                * d1**2
                * kj**3
                * (1 + K)
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                - 16
                * d1**2
                * kj**4
                * (-1 + 2 * K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                - 2
                * kj**2
                * (
                    2 * d1**2 * (5 + 2 * K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (4 + 3 * K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    - 6
                    * K
                    * (3 + K)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + P022
                * (
                    (4 * d1**2 * kj**2 * (5 + 2 * K) + 3 * (24 + 11 * K + K**2))
                    * P041
                    * Pf1
                    - 3 * (24 + 11 * K + K**2) * P031 * Pf2
                    - 2
                    * d1
                    * kj**2
                    * (
                        -3 * (4 + 3 * K) * P043 * Pf1
                        + (2 * d1 * (5 + 2 * K) * P031 + 3 * (4 + 3 * K) * P033) * Pf2
                    )
                    - 16 * d1**2 * kj**3 * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                    - 6
                    * kj
                    * (3 + K)
                    * (
                        (3 + K) * P043 * Pf1
                        + P044 * Pf1
                        - ((3 + K) * P033 + P034) * Pf2
                    )
                )
                + 6
                * kj
                * (3 + K)
                * (
                    P024
                    * (
                        (-(3 + K)) * P041 * Pf1
                        + P042 * Pf1
                        + ((3 + K) * P031 - P032) * Pf2
                    )
                    + P023
                    * (
                        (-P041) * Pf1
                        + (3 + K) * P042 * Pf1
                        + (P031 - (3 + K) * P032) * Pf2
                    )
                    + P021
                    * (
                        P043 * Pf1
                        + (3 + K) * P044 * Pf1
                        - (P033 + (3 + K) * P034) * Pf2
                    )
                )
            )
            + C1
            * (
                (108 + 57 * K + 7 * K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                - 4
                * kj**2
                * (18 + 21 * K + 5 * K**2)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + 2
                * kj
                * (
                    P023
                    * (
                        (18 + 9 * K + K**2) * P041 * Pf1
                        - 6 * (5 + 7 * K + K**2) * P042 * Pf1
                        - ((18 + 9 * K + K**2) * P031 - 6 * (5 + 7 * K + K**2) * P032)
                        * Pf2
                    )
                    + P024
                    * (
                        6 * (5 + 7 * K + K**2) * P041 * Pf1
                        - (18 + 9 * K + K**2) * P042 * Pf1
                        - (6 * (5 + 7 * K + K**2) * P031 - (18 + 9 * K + K**2) * P032)
                        * Pf2
                    )
                    + P021
                    * (
                        (-(18 + 9 * K + K**2)) * P043 * Pf1
                        - 6 * (5 + 7 * K + K**2) * P044 * Pf1
                        + ((18 + 9 * K + K**2) * P033 + 6 * (5 + 7 * K + K**2) * P034)
                        * Pf2
                    )
                )
                + P022
                * (
                    (-(108 + 57 * K + 7 * K**2)) * P041 * Pf1
                    + (108 + 57 * K + 7 * K**2) * P031 * Pf2
                    + 2
                    * kj
                    * (
                        6 * (5 + 7 * K + K**2) * P043 * Pf1
                        + (18 + 9 * K + K**2) * P044 * Pf1
                        - (6 * (5 + 7 * K + K**2) * P033 + (18 + 9 * K + K**2) * P034)
                        * Pf2
                    )
                )
            )
            * S1
            + 2
            * (
                -2 * (1 + 8 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                + 8
                * d1**2
                * kj**3
                * (1 + K)
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                + 8
                * d1**2
                * kj**4
                * (-1 + 2 * K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + kj**2
                * (
                    2 * d1**2 * (5 + 2 * K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (4 + 3 * K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + 4
                    * (1 + 8 * K + K**2)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + kj
                * (
                    P023
                    * (
                        (-(1 + 8 * K + K**2)) * P041 * Pf1
                        + 3 * (3 + K) ** 2 * P042 * Pf1
                        + ((1 + 8 * K + K**2) * P031 - 3 * (3 + K) ** 2 * P032) * Pf2
                    )
                    + P024
                    * (
                        -3 * (3 + K) ** 2 * P041 * Pf1
                        + (1 + 8 * K + K**2) * P042 * Pf1
                        + (3 * (3 + K) ** 2 * P031 - (1 + 8 * K + K**2) * P032) * Pf2
                    )
                    + P021
                    * (
                        (1 + 8 * K + K**2) * P043 * Pf1
                        + 3 * (3 + K) ** 2 * P044 * Pf1
                        - ((1 + 8 * K + K**2) * P033 + 3 * (3 + K) ** 2 * P034) * Pf2
                    )
                )
                + P022
                * (
                    2 * (1 + 8 * K + K**2 - d1**2 * kj**2 * (5 + 2 * K)) * P041 * Pf1
                    - 2 * (1 + 8 * K + K**2) * P031 * Pf2
                    + d1
                    * kj**2
                    * (
                        -3 * (4 + 3 * K) * P043 * Pf1
                        + (2 * d1 * (5 + 2 * K) * P031 + 3 * (4 + 3 * K) * P033) * Pf2
                    )
                    + 8 * d1**2 * kj**3 * (1 + K) * (P044 * Pf1 - P034 * Pf2)
                    + kj
                    * (
                        -3 * (3 + K) ** 2 * P043 * Pf1
                        - (1 + 8 * K + K**2) * P044 * Pf1
                        + (3 * (3 + K) ** 2 * P033 + (1 + 8 * K + K**2) * P034) * Pf2
                    )
                )
            )
            * S1**2
        )
        + B1**2
        * (
            C1**2
            * (
                -2 * (9 + 12 * K + K**2) * P021 * (P042 * Pf1 - P032 * Pf2)
                + 8
                * d1**2
                * kj**3
                * (4 + 3 * K)
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                + 16
                * d1**2
                * kj**4
                * (1 + 3 * K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                - kj
                * (9 + 12 * K + K**2)
                * (
                    P024 * (3 * P041 * Pf1 - P042 * Pf1 - 3 * P031 * Pf2 + P032 * Pf2)
                    + P023 * (P041 * Pf1 - 3 * P042 * Pf1 - P031 * Pf2 + 3 * P032 * Pf2)
                    + P021
                    * ((-P043) * Pf1 - 3 * P044 * Pf1 + P033 * Pf2 + 3 * P034 * Pf2)
                )
                + 2
                * kj**2
                * (
                    2 * d1**2 * (7 + 3 * K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (8 + 5 * K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + 2
                    * (9 + 12 * K + K**2)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + P022
                * (
                    2
                    * (9 + 12 * K + K**2 - 2 * d1**2 * kj**2 * (7 + 3 * K))
                    * P041
                    * Pf1
                    - 2 * (9 + 12 * K + K**2) * P031 * Pf2
                    + 2
                    * d1
                    * kj**2
                    * (
                        -3 * (8 + 5 * K) * P043 * Pf1
                        + (2 * d1 * (7 + 3 * K) * P031 + 3 * (8 + 5 * K) * P033) * Pf2
                    )
                    + 8 * d1**2 * kj**3 * (4 + 3 * K) * (P044 * Pf1 - P034 * Pf2)
                    - kj
                    * (9 + 12 * K + K**2)
                    * (3 * P043 * Pf1 + P044 * Pf1 - (3 * P033 + P034) * Pf2)
                )
            )
            - 3
            * C1
            * (
                (-(36 + 21 * K + K**2)) * P021 * (P042 * Pf1 - P032 * Pf2)
                + 4
                * kj**2
                * K
                * (5 + K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + P022
                * (
                    (36 + 21 * K + K**2) * P041 * Pf1
                    - (36 + 21 * K + K**2) * P031 * Pf2
                    - 2
                    * kj
                    * (
                        (13 + 11 * K + K**2) * P043 * Pf1
                        + (3 + 2 * K) * P044 * Pf1
                        - ((13 + 11 * K + K**2) * P033 + (3 + 2 * K) * P034) * Pf2
                    )
                )
                - 2
                * kj
                * (
                    P024
                    * (
                        (13 + 11 * K + K**2) * P041 * Pf1
                        - (3 + 2 * K) * P042 * Pf1
                        - ((13 + 11 * K + K**2) * P031 - (3 + 2 * K) * P032) * Pf2
                    )
                    + P023
                    * (
                        (3 + 2 * K) * P041 * Pf1
                        - (13 + 11 * K + K**2) * P042 * Pf1
                        + ((-(3 + 2 * K)) * P031 + (13 + 11 * K + K**2) * P032) * Pf2
                    )
                    + P021
                    * (
                        (-(3 + 2 * K)) * P043 * Pf1
                        - (13 + 11 * K + K**2) * P044 * Pf1
                        + ((3 + 2 * K) * P033 + (13 + 11 * K + K**2) * P034) * Pf2
                    )
                )
            )
            * S1
            + (
                (-(106 + 35 * K + K**2)) * P021 * (P042 * Pf1 - P032 * Pf2)
                - 8
                * d1**2
                * kj**3
                * (4 + 3 * K)
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                - 16
                * d1**2
                * kj**4
                * (1 + 3 * K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                - 2
                * kj**2
                * (
                    2 * d1**2 * (7 + 3 * K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (8 + 5 * K)
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    - 2
                    * (-1 + K + 2 * K**2)
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + P022
                * (
                    (106 + 35 * K + K**2 + 4 * d1**2 * kj**2 * (7 + 3 * K)) * P041 * Pf1
                    - (106 + 35 * K + K**2) * P031 * Pf2
                    - 2
                    * d1
                    * kj**2
                    * (
                        -3 * (8 + 5 * K) * P043 * Pf1
                        + (2 * d1 * (7 + 3 * K) * P031 + 3 * (8 + 5 * K) * P033) * Pf2
                    )
                    - 8 * d1**2 * kj**3 * (4 + 3 * K) * (P044 * Pf1 - P034 * Pf2)
                    - kj
                    * (
                        3 * (9 + 12 * K + K**2) * P043 * Pf1
                        - (-17 + 2 * K + K**2) * P044 * Pf1
                        - (3 * (9 + 12 * K + K**2) * P033 - (-17 + 2 * K + K**2) * P034)
                        * Pf2
                    )
                )
                + kj
                * (
                    P024
                    * (
                        -3 * (9 + 12 * K + K**2) * P041 * Pf1
                        - (-17 + 2 * K + K**2) * P042 * Pf1
                        + (3 * (9 + 12 * K + K**2) * P031 + (-17 + 2 * K + K**2) * P032)
                        * Pf2
                    )
                    + P023
                    * (
                        (-17 + 2 * K + K**2) * P041 * Pf1
                        + 3 * (9 + 12 * K + K**2) * P042 * Pf1
                        - ((-17 + 2 * K + K**2) * P031 + 3 * (9 + 12 * K + K**2) * P032)
                        * Pf2
                    )
                    + P021
                    * (
                        (-(-17 + 2 * K + K**2)) * P043 * Pf1
                        + 3 * (9 + 12 * K + K**2) * P044 * Pf1
                        + ((-17 + 2 * K + K**2) * P033 - 3 * (9 + 12 * K + K**2) * P034)
                        * Pf2
                    )
                )
            )
            * S1**2
        )
        + B2**2
        * (
            C1**2
            * (
                -12
                * d1
                * kj**2
                * (2 + K)
                * (
                    P024 * P041 * Pf1
                    + P023 * P042 * Pf1
                    - P021 * P044 * Pf1
                    - P024 * P031 * Pf2
                    - P023 * P032 * Pf2
                    + P021 * P034 * Pf2
                )
                - 9
                * (
                    (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + kj
                    * K
                    * (
                        P024 * P041 * Pf1
                        - P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        + P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                )
                + 4
                * d1**2
                * kj**2
                * (2 + K)
                * (
                    (-P021) * P042 * Pf1
                    + P021 * P032 * Pf2
                    - 2
                    * kj
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    - 4
                    * kj**2
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + P022
                * (
                    (9 + 4 * d1**2 * kj**2) * (2 + K) * P041 * Pf1
                    + 12 * d1 * kj**2 * (2 + K) * (P043 * Pf1 - P033 * Pf2)
                    - 4
                    * d1**2
                    * kj**2
                    * (2 + K)
                    * (2 * kj * P044 * Pf1 + P031 * Pf2 - 2 * kj * P034 * Pf2)
                    - 9 * ((2 + K) * P031 * Pf2 + kj * K * (P043 * Pf1 - P033 * Pf2))
                )
            )
            - 6
            * C1
            * (
                -4
                * kj
                * (
                    P024 * P041 * Pf1
                    - P023 * P042 * Pf1
                    + P022 * P043 * Pf1
                    - P021 * P044 * Pf1
                    - P024 * P031 * Pf2
                    + P023 * P032 * Pf2
                    - P022 * P033 * Pf2
                    + P021 * P034 * Pf2
                )
                + K
                * (
                    -2 * P021 * P042 * Pf1
                    + 2 * P021 * P032 * Pf2
                    + 4
                    * kj**2
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + P022
                    * (
                        2 * P041 * Pf1
                        - 2 * P031 * Pf2
                        + kj
                        * (-2 * P043 * Pf1 - P044 * Pf1 + 2 * P033 * Pf2 + P034 * Pf2)
                    )
                    + kj
                    * (
                        P023
                        * ((-P041) * Pf1 + 2 * P042 * Pf1 + (P031 - 2 * P032) * Pf2)
                        + P024
                        * (-2 * P041 * Pf1 + P042 * Pf1 + 2 * P031 * Pf2 - P032 * Pf2)
                        + P021 * (P043 * Pf1 + 2 * P044 * Pf1 - (P033 + 2 * P034) * Pf2)
                    )
                )
            )
            * S1
            + (
                (2 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                + 8
                * d1**2
                * kj**3
                * (2 + K)
                * (
                    P023 * P041 * Pf1
                    - P024 * P042 * Pf1
                    - P021 * P043 * Pf1
                    - P023 * P031 * Pf2
                    + P024 * P032 * Pf2
                    + P021 * P033 * Pf2
                )
                + 16
                * d1**2
                * kj**4
                * (2 + K)
                * (
                    P024 * P043 * Pf1
                    - P023 * P044 * Pf1
                    - P024 * P033 * Pf2
                    + P023 * P034 * Pf2
                )
                + P022
                * (
                    (-(1 + 4 * d1**2 * kj**2)) * (2 + K) * P041 * Pf1
                    - 9 * kj * K * P043 * Pf1
                    - 8 * kj * P044 * Pf1
                    - 4 * kj * K * P044 * Pf1
                    + 2 * P031 * Pf2
                    + K * P031 * Pf2
                    + 9 * kj * K * P033 * Pf2
                    + 8 * kj * P034 * Pf2
                    + 4 * kj * K * P034 * Pf2
                    - 12 * d1 * kj**2 * (2 + K) * (P043 * Pf1 - P033 * Pf2)
                    + 4
                    * d1**2
                    * kj**2
                    * (2 + K)
                    * (2 * kj * P044 * Pf1 + P031 * Pf2 - 2 * kj * P034 * Pf2)
                )
                + 4
                * kj**2
                * (2 + K)
                * (
                    d1**2 * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 3
                    * d1
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + 4
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                )
                + kj
                * (
                    P023
                    * (
                        -4 * (2 + K) * P041 * Pf1
                        + 9 * K * P042 * Pf1
                        + 8 * P031 * Pf2
                        + 4 * K * P031 * Pf2
                        - 9 * K * P032 * Pf2
                    )
                    + 8
                    * (
                        P024 * P042 * Pf1
                        + P021 * P043 * Pf1
                        - P024 * P032 * Pf2
                        - P021 * P033 * Pf2
                    )
                    + K
                    * (
                        P024
                        * (
                            -9 * P041 * Pf1
                            + 4 * P042 * Pf1
                            + 9 * P031 * Pf2
                            - 4 * P032 * Pf2
                        )
                        + P021
                        * (
                            4 * P043 * Pf1
                            + 9 * P044 * Pf1
                            - 4 * P033 * Pf2
                            - 9 * P034 * Pf2
                        )
                    )
                )
            )
            * S1**2
        )
    )

    N7 = (
        -3
        * (C1 - S1)
        * (
            2
            * B1
            * (
                2
                * C1
                * (
                    2 * (3 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 4
                    * d1**2
                    * kj**3
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    - 8
                    * d1**2
                    * kj**4
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + kj
                    * (3 + K)
                    * (
                        P024
                        * (3 * P041 * Pf1 - P042 * Pf1 - 3 * P031 * Pf2 + P032 * Pf2)
                        + P023
                        * (P041 * Pf1 - 3 * P042 * Pf1 - P031 * Pf2 + 3 * P032 * Pf2)
                        + P021
                        * ((-P043) * Pf1 - 3 * P044 * Pf1 + P033 * Pf2 + 3 * P034 * Pf2)
                    )
                    - 2
                    * kj**2
                    * (
                        d1**2 * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 3
                        * d1
                        * (
                            P024 * P041 * Pf1
                            + P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            - P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                        + 2
                        * (3 + K)
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        2 * (-3 + d1**2 * kj**2 - K) * P041 * Pf1
                        + 2 * (3 + K) * P031 * Pf2
                        - 2
                        * d1
                        * kj**2
                        * (-3 * P043 * Pf1 + d1 * P031 * Pf2 + 3 * P033 * Pf2)
                        + 4 * d1**2 * kj**3 * ((-P044) * Pf1 + P034 * Pf2)
                        + kj
                        * (3 + K)
                        * (3 * P043 * Pf1 + P044 * Pf1 - (3 * P033 + P034) * Pf2)
                    )
                )
                + (
                    -3 * (8 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    - 8
                    * d1**2
                    * kj**3
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    - 16
                    * d1**2
                    * kj**4
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    - 4
                    * kj**2
                    * (
                        d1**2 * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 3
                        * d1
                        * (
                            P024 * P041 * Pf1
                            + P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            - P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                        + 3
                        * K
                        * (
                            (-P024) * P043 * Pf1
                            + P023 * P044 * Pf1
                            + P024 * P033 * Pf2
                            - P023 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        (4 * d1**2 * kj**2 + 3 * (8 + K)) * P041 * Pf1
                        - 3 * (8 + K) * P031 * Pf2
                        - 4
                        * d1
                        * kj**2
                        * (-3 * P043 * Pf1 + d1 * P031 * Pf2 + 3 * P033 * Pf2)
                        - 8 * d1**2 * kj**3 * (P044 * Pf1 - P034 * Pf2)
                        - 6
                        * kj
                        * (
                            (3 + K) * P043 * Pf1
                            + P044 * Pf1
                            - ((3 + K) * P033 + P034) * Pf2
                        )
                    )
                    - 6
                    * kj
                    * (
                        P024
                        * (
                            (3 + K) * P041 * Pf1
                            - P042 * Pf1
                            - 3 * P031 * Pf2
                            - K * P031 * Pf2
                            + P032 * Pf2
                        )
                        + P023
                        * (
                            P041 * Pf1
                            - (3 + K) * P042 * Pf1
                            + (-P031 + (3 + K) * P032) * Pf2
                        )
                        + P021
                        * (
                            (-P043) * Pf1
                            - (3 + K) * P044 * Pf1
                            + (P033 + (3 + K) * P034) * Pf2
                        )
                    )
                )
                * S1
            )
            + B2
            * (
                C1
                * (
                    3 * (8 + K) * P021 * (P042 * Pf1 - P032 * Pf2)
                    + 16
                    * d1**2
                    * kj**3
                    * (
                        P023 * P041 * Pf1
                        - P024 * P042 * Pf1
                        - P021 * P043 * Pf1
                        - P023 * P031 * Pf2
                        + P024 * P032 * Pf2
                        + P021 * P033 * Pf2
                    )
                    + 32
                    * d1**2
                    * kj**4
                    * (
                        P024 * P043 * Pf1
                        - P023 * P044 * Pf1
                        - P024 * P033 * Pf2
                        + P023 * P034 * Pf2
                    )
                    + 4
                    * kj**2
                    * (
                        2 * d1**2 * P021 * (P042 * Pf1 - P032 * Pf2)
                        + 6
                        * d1
                        * (
                            P024 * P041 * Pf1
                            + P023 * P042 * Pf1
                            - P021 * P044 * Pf1
                            - P024 * P031 * Pf2
                            - P023 * P032 * Pf2
                            + P021 * P034 * Pf2
                        )
                        + 3
                        * K
                        * (
                            (-P024) * P043 * Pf1
                            + P023 * P044 * Pf1
                            + P024 * P033 * Pf2
                            - P023 * P034 * Pf2
                        )
                    )
                    + P022
                    * (
                        (-(8 * d1**2 * kj**2 + 3 * (8 + K))) * P041 * Pf1
                        + 3 * (8 + K) * P031 * Pf2
                        + 8
                        * d1
                        * kj**2
                        * (-3 * P043 * Pf1 + d1 * P031 * Pf2 + 3 * P033 * Pf2)
                        + 16 * d1**2 * kj**3 * (P044 * Pf1 - P034 * Pf2)
                        + 6
                        * kj
                        * (
                            (3 + K) * P043 * Pf1
                            + P044 * Pf1
                            - ((3 + K) * P033 + P034) * Pf2
                        )
                    )
                    + 6
                    * kj
                    * (
                        P024
                        * (
                            (3 + K) * P041 * Pf1
                            - P042 * Pf1
                            - 3 * P031 * Pf2
                            - K * P031 * Pf2
                            + P032 * Pf2
                        )
                        + P023
                        * (
                            P041 * Pf1
                            - (3 + K) * P042 * Pf1
                            + (-P031 + (3 + K) * P032) * Pf2
                        )
                        + P021
                        * (
                            (-P043) * Pf1
                            - (3 + K) * P044 * Pf1
                            + (P033 + (3 + K) * P034) * Pf2
                        )
                    )
                )
                + (
                    24
                    * d1
                    * kj**2
                    * (
                        P024 * P041 * Pf1
                        + P023 * P042 * Pf1
                        - P022 * P043 * Pf1
                        - P021 * P044 * Pf1
                        - P024 * P031 * Pf2
                        - P023 * P032 * Pf2
                        + P022 * P033 * Pf2
                        + P021 * P034 * Pf2
                    )
                    + 6
                    * kj
                    * (
                        P021 * P043 * Pf1
                        - 3 * P022 * P043 * Pf1
                        + 3 * P021 * P044 * Pf1
                        - P022 * P044 * Pf1
                        - P021 * P033 * Pf2
                        + 3 * P022 * P033 * Pf2
                        - 3 * P021 * P034 * Pf2
                        + P022 * P034 * Pf2
                        + P024
                        * (
                            -3 * P041 * Pf1
                            + P042 * Pf1
                            + 8 * kj * P043 * Pf1
                            + 3 * P031 * Pf2
                            - P032 * Pf2
                            - 8 * kj * P033 * Pf2
                        )
                        + P023
                        * (
                            (-P041) * Pf1
                            + 3 * P042 * Pf1
                            - 8 * kj * P044 * Pf1
                            + P031 * Pf2
                            - 3 * P032 * Pf2
                            + 8 * kj * P034 * Pf2
                        )
                    )
                    - 8
                    * d1**2
                    * kj**2
                    * (
                        (-P021) * P042 * Pf1
                        + P021 * P032 * Pf2
                        - 2
                        * kj
                        * (
                            P023 * P041 * Pf1
                            - P024 * P042 * Pf1
                            - P021 * P043 * Pf1
                            - P023 * P031 * Pf2
                            + P024 * P032 * Pf2
                            + P021 * P033 * Pf2
                        )
                        + P022
                        * (
                            P041 * Pf1
                            - 2 * kj * P044 * Pf1
                            - P031 * Pf2
                            + 2 * kj * P034 * Pf2
                        )
                        - 4
                        * kj**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                    )
                    + K
                    * (
                        -5 * P021 * P042 * Pf1
                        + 5 * P021 * P032 * Pf2
                        + 4
                        * kj**2
                        * (
                            P024 * P043 * Pf1
                            - P023 * P044 * Pf1
                            - P024 * P033 * Pf2
                            + P023 * P034 * Pf2
                        )
                        + kj
                        * (
                            P023
                            * (
                                -4 * P041 * Pf1
                                + 6 * P042 * Pf1
                                + 4 * P031 * Pf2
                                - 6 * P032 * Pf2
                            )
                            + P024
                            * (
                                -6 * P041 * Pf1
                                + 4 * P042 * Pf1
                                + 6 * P031 * Pf2
                                - 4 * P032 * Pf2
                            )
                            + 2
                            * P021
                            * (
                                2 * P043 * Pf1
                                + 3 * P044 * Pf1
                                - 2 * P033 * Pf2
                                - 3 * P034 * Pf2
                            )
                        )
                        + P022
                        * (
                            5 * P041 * Pf1
                            - 5 * P031 * Pf2
                            + kj
                            * (
                                -6 * P043 * Pf1
                                - 4 * P044 * Pf1
                                + 6 * P033 * Pf2
                                + 4 * P034 * Pf2
                            )
                        )
                    )
                )
                * S1
            )
        )
    )

    N8 = (
        9
        * (
            -2 * P021 * P042 * Pf1
            + 2 * P021 * P032 * Pf2
            + 4
            * kj**2
            * (
                P024 * P043 * Pf1
                - P023 * P044 * Pf1
                - P024 * P033 * Pf2
                + P023 * P034 * Pf2
            )
            + P022
            * (
                2 * P041 * Pf1
                - 2 * P031 * Pf2
                + kj * (-3 * P043 * Pf1 - P044 * Pf1 + 3 * P033 * Pf2 + P034 * Pf2)
            )
            + kj
            * (
                P023 * ((-P041) * Pf1 + 3 * P042 * Pf1 + (P031 - 3 * P032) * Pf2)
                + P024 * (-3 * P041 * Pf1 + P042 * Pf1 + 3 * P031 * Pf2 - P032 * Pf2)
                + P021 * (P043 * Pf1 + 3 * P044 * Pf1 - (P033 + 3 * P034) * Pf2)
            )
        )
        * (C1 - S1) ** 2
    )

    return [N1, N2, N3, N4, N5, N6, N7, N8]


def roots_2visc(P0, Pzsf, d1, kj, B1, B2, K):
    N1 = get_num1_2visc(P0, Pzsf, d1, kj, B1, B2, K)
    N2 = get_num2_2visc(P0, Pzsf, d1, kj, B1, B2, K)
    [D1, D2, D3, D4, D5, D6, D7, D8] = get_denom_2visc(P0, d1, kj, B1, B2, K)

    # insert a zero at the beginning of D to set up the polynomial
    D = [0, D1, D2, D3, D4, D5, D6, D7, D8]
    roots = np.polynomial.Polynomial(D).roots()
    return N1, N2, roots, D8


# for 0 order up to +/- 2nd order, feed in point sources, propagator matrix, Num1/2 matrix, and viscous roots
# to obtain displacements
def get_UsUn4(F4, Pzs4, nu, roots, a, T, t_eq):

    nu1, nu2 = np.dot(Pzs4[2:4, :], F4) @ nu

    US = []
    UR = []
    for idx in range(len(roots) - 1):
        s = roots[idx]
        mask = np.ones_like(roots, dtype=bool)
        mask[[idx, -1]] = False
        Den = a * s * np.prod(s - roots[mask])
        Num1 = -np.dot(nu1, s ** np.arange(len(nu1)))
        Num2 = -np.dot(nu2, s ** np.arange(len(nu2)))
        Tterm = 1 / (1 - np.exp(s * T))
        term1 = T * Num1 * s * np.exp(s * t_eq) * Tterm / Den + Num1 / Den
        term2 = T * Num2 * s * np.exp(s * t_eq) * Tterm / Den + Num2 / Den

        US.append(term1)
        UR.append(term2)

    US = np.array(US).sum()
    UR = -np.array(UR).sum()

    return US, UR


def antiplane_2visc_terms(PzsF, B1, B2, k, d1, P0, mu):
    Ph11, Ph12 = P0[0, 0], P0[0, 1]
    Ph21, Ph22 = P0[1, 0], P0[1, 1]

    C1 = np.cosh(abs(k) * d1)
    Z1 = np.sinh(abs(k) * d1)

    Pf2 = PzsF[1]

    N1 = -4 * B1 * B2 * C1 * Pf2 * Ph11 + 4 * B1**2 * Pf2 * Ph11 * Z1

    N2 = (
        -2 * B1 * C1 * Pf2 * Ph11
        - 2 * B2 * C1 * Pf2 * Ph11
        + 4 * B1 * Pf2 * Ph11 * Z1
        + 2 * B1 * C1 * Pf2 * Ph12 * mu * abs(k)
        - 2 * B2 * Pf2 * Ph12 * mu * Z1 * abs(k)
    )

    N3 = (
        -C1 * Pf2 * Ph11
        + Pf2 * Ph11 * Z1
        + C1 * Pf2 * Ph12 * mu * abs(k)
        - Pf2 * Ph12 * mu * Z1 * abs(k)
    )

    D1 = 4 * B1 * B2 * C1 * Ph21 - 4 * B1**2 * Ph21 * Z1

    D2 = (
        2 * B1 * C1 * Ph21
        + 2 * B2 * C1 * Ph21
        - 4 * B1 * Ph21 * Z1
        - 2 * B1 * C1 * Ph22 * mu * abs(k)
        + 2 * B2 * Ph22 * mu * Z1 * abs(k)
    )

    D3 = C1 * Ph21 - Ph21 * Z1 - C1 * Ph22 * mu * abs(k) + Ph22 * mu * Z1 * abs(k)

    return N1, N2, N3, D1, D2, D3


def roots_2visc_antiplane(PzsF, B1, B2, k, d1, P0, u):
    N1, N2, N3, D1, D2, D3 = antiplane_2visc_terms(PzsF, B1, B2, k, d1, P0, u)

    R = np.polynomial.Polynomial([0, D1, D2, D3]).roots()

    return N1, N2, N3, R, D3


def get_2d_U(F0, Prop, N, R, a, T, t):
    N1, N2, N3 = np.dot(N, (Prop @ F0))

    s = R[0]
    Den = a * s * (s - R[1])
    Num = -(N1 + N2 * s + N3 * s**2)
    Tterm = 1 / (1 - np.exp(s * T))
    term1 = T * Num * s * np.exp(s * t) * Tterm / Den + Num / Den

    s = R[1]
    Den = a * s * (s - R[0])
    Num = -(N1 + N2 * s + N3 * s**2)
    Tterm = 1 / (1 - np.exp(s * T))
    term2 = T * Num * s * np.exp(s * t) * Tterm / Den + Num / Den

    UT = -(term1 + term2)

    return UT


# FIX this. we loop through it a bunch of times, just gfeed it constants
# instead of running get_constants() hundreds of times
def get_bessel_coeffs(
    Ph2, Ph4, Pzs4, Pzs2, M, k, lam, mu, beta1, beta2, H1, H2, timeEq, T
):

    consts = get_constants()

    Mxx, Myy, Mzz = M[0, 0], M[1, 1], M[2, 2]
    Mxy, Mxz, Myz = M[0, 1], M[0, 2], M[1, 2]

    [Nu11, Nu21, R4x4, a4x4] = roots_2visc(Ph4, [1, 0], consts['dH'], k, beta1, beta2, consts['bulk'])
    [Nu12, Nu22, R4x4, a4x4] = roots_2visc(Ph4, [0, 1], consts['dH'], k, beta1, beta2, consts['bulk'])

    Nu = np.array([[Nu11, Nu21], [Nu12, Nu22]])

    [N11, N21, N31, R2x2, a2x2] = roots_2visc_antiplane(
        [1, 0], beta1, beta2, k, consts['dH'], Ph2, consts['Gshear']
    )
    [N12, N22, N32, R2x2, a2x2] = roots_2visc_antiplane(
        [0, 1], beta1, beta2, k, consts["dH"], Ph2, consts["Gshear"]
    )

    N_2x2 = np.array([[N11, N12], [N21, N22], [N31, N32]])

    # Order = 0
    F4x4_0 = np.array([0, Mzz / (lam + 2 * mu),
            k * ((-Mxx - Myy) / 2 + lam * Mzz / (lam + 2 * mu)),0,])

    # Order = -1
    F4x4_n1 = np.array([-1 / 2 * (Mxz + 1j * Myz) / mu, 0, 0, 0])
    F2x2_n1 = np.array([1 / 2 * (Myz - 1j * Mxz) / mu, 0])

    # Order = 1
    F4x4_p1 = np.array([-1 / 2 * (-Mxz + 1j * Myz) / mu, 0, 0, 0])
    F2x2_p1 = np.array([-1 / 2 * (Myz + 1j * Mxz) / mu, 0])

    # Order = -2
    F4x4_n2 = np.array([0, 0, k * (-(Myy - Mxx) / 4 + 1j * Mxy / 2), 0])
    F2x2_n2 = np.array([0, k * (1j * (Mxx - Myy) / 4 - 1 / 2 * Mxy)])

    # Order = +2
    F4x4_p2 = np.array([0, 0, k * (-(Myy - Mxx) / 4 - 1j * Mxy / 2), 0])
    F2x2_p2 = np.array([0, k * (-1j * (Mxx - Myy) / 4 - 1 / 2 * Mxy)])

    # for 0 order up to +/- 2nd order, feed in point sources, propagator matrix, Num1/2 matrix, and viscous roots
    US_o, UR_o = get_UsUn4(F4x4_0, Pzs4, Nu, R4x4, a4x4, T, timeEq)
    US_n1, UR_n1 = get_UsUn4(F4x4_n1, Pzs4, Nu, R4x4, a4x4, T, timeEq)
    US_p1, UR_p1 = get_UsUn4(F4x4_p1, Pzs4, Nu, R4x4, a4x4, T, timeEq)
    US_n2, UR_n2 = get_UsUn4(F4x4_n2, Pzs4, Nu, R4x4, a4x4, T, timeEq)
    US_p2, UR_p2 = get_UsUn4(F4x4_p2, Pzs4, Nu, R4x4, a4x4, T, timeEq)

    # for orders +/-1 and 2, get 2x2 solution too
    UT_n1 = get_2d_U(F2x2_n1, Pzs2, N_2x2, R2x2, a2x2, T, timeEq)
    UT_p1 = get_2d_U(F2x2_p1, Pzs2, N_2x2, R2x2, a2x2, T, timeEq)
    UT_n2 = get_2d_U(F2x2_n2, Pzs2, N_2x2, R2x2, a2x2, T, timeEq)
    UT_p2 = get_2d_U(F2x2_p2, Pzs2, N_2x2, R2x2, a2x2, T, timeEq)

    return (
        UR_o,
        US_o,
        UR_n1,
        US_n1,
        UT_n1,
        UR_p1,
        US_p1,
        UT_p1,
        UR_n2,
        US_n2,
        UT_n2,
        UR_p2,
        US_p2,
        UT_p2,
    )


def bessel_wrapper(zpos, D, mu, lam, M, beta1, beta2, H1, H2, timeEq, T, k):

    Us_i = []
    for z_o in zpos:

        slip_vec = []

        for k_o in k:

            Pzs4, Pzs2, Ph4, Ph2 = get_prop(H1, k_o, mu, lam, D + z_o)

            slip_vec.append(
                get_bessel_coeffs(
                    Ph2,
                    Ph4,
                    Pzs4,
                    Pzs2,
                    M,
                    k_o,
                    lam,
                    mu,
                    beta1,
                    beta2,
                    H1,
                    H2,
                    timeEq,
                    T,
                )
            )

        Us_i.append(slip_vec)

    Us_i = np.array(Us_i)

    return Us_i


def integrate_disps(U_comp, k, k_tile, r_tile, theta, b0, b1):

    b2 = -b0 + 2 / (k_tile * r_tile) * b1
    b3 = -b1 + 4 / (k_tile * r_tile) * b2

    theta_tile = np.tile(theta, (len(k), 1)).T

    nobs = k_tile.shape[0]

    # break out into individual compoennts for legibility
    UR_o0 = np.tile(U_comp[:, 0], (nobs, 1))
    US_o0 = np.tile(U_comp[:, 1], (nobs, 1))

    UR_n1 = np.tile(U_comp[:, 2], (nobs, 1))
    US_n1 = np.tile(U_comp[:, 3], (nobs, 1))
    UT_n1 = np.tile(U_comp[:, 4], (nobs, 1))

    UR_p1 = np.tile(U_comp[:, 5], (nobs, 1))
    US_p1 = np.tile(U_comp[:, 6], (nobs, 1))
    UT_p1 = np.tile(U_comp[:, 7], (nobs, 1))

    UR_n2 = np.tile(U_comp[:, 8], (nobs, 1))
    US_n2 = np.tile(U_comp[:, 9], (nobs, 1))
    UT_n2 = np.tile(U_comp[:, 10], (nobs, 1))

    UR_p2 = np.tile(U_comp[:, 11], (nobs, 1))
    US_p2 = np.tile(U_comp[:, 12], (nobs, 1))
    UT_p2 = np.tile(U_comp[:, 13], (nobs, 1))

    # preallocate slip components
    Uz = np.ones((nobs, 5))
    Ur = np.ones((nobs, 5))
    Utheta = np.zeros((nobs, 5))

    # calculate 0th order components...
    M = 0
    DrJm = -k_tile * b1

    Uz_integrand = k_tile * UR_o0 * b0
    Ur_integrand = US_o0 * DrJm

    Uz[:, 0] = np.real((1 / (2 * pi)) * np.trapz(Uz_integrand, k, 2))
    Ur[:, 0] = np.real((1 / (2 * pi)) * np.trapz(Ur_integrand, k, 2))
    # Utheta2 is 0 for 0th order

    #  calculate -1st order components...
    M = -1
    DrJm = 0.5 * k_tile * (b2 - b0)

    Uz_integrand = -k_tile * UR_n1 * b1 * np.exp(1j * M * theta_tile)
    Ur_integrand = (US_n1 * DrJm - (0 + 1j) * M * UT_n1 * (1 / r_tile) * b1) * np.exp(
        1j * M * theta_tile
    )

    Uz[:, 1] = np.real((1 / (2 * pi)) * np.trapz(Uz_integrand, k, 2))
    Ur[:, 1] = np.real((1 / (2 * pi)) * np.trapz(Ur_integrand, k, 2))
    # Utheta is 0 for -1st order

    # calculate 1st order components...
    M = 1
    DrJm = 0.5 * k_tile * (b0 - b2)

    Uz_integrand = k_tile * UR_p1 * b1 * np.exp(1j * M * theta_tile)
    Ur_integrand = (US_p1 * DrJm + (0 + 1j) * M * UT_p1 * (1 / r_tile) * b1) * np.exp(
        1j * M * theta_tile
    )

    Uz[:, 2] = np.real((1 / (2 * pi)) * np.trapz(Uz_integrand, k, 2))
    Ur[:, 2] = np.real((1 / (2 * pi)) * np.trapz(Ur_integrand, k, 2))
    # Utheta is 0 for +1st order

    # calculate -2nd order components...
    M = -2
    DrJm = 0.5 * k_tile * (b1 - b3)

    Uz_integrand = k_tile * UR_n2 * b2 * np.exp(1j * M * theta_tile)
    Ur_integrand = (US_n2 * DrJm + (0 + 1j) * M * UT_n2 * (1 / r_tile) * b2) * np.exp(
        1j * M * theta_tile
    )
    Utheta_integrand = (1j * M * US_n2 * (1 / r_tile) * b2 - UT_n2 * DrJm) * np.exp(
        1j * M * theta_tile
    )

    Uz[:, 3] = np.real((1 / (2 * pi)) * np.trapz(Uz_integrand, k, 2))
    Ur[:, 3] = np.real((1 / (2 * pi)) * np.trapz(Ur_integrand, k, 2))
    Utheta[:, 3] = np.real((1 / (2 * pi)) * np.trapz(Utheta_integrand, k, 2))

    # calculate 2nd order components...
    M = 2
    DrJm = 0.5 * k_tile * (b1 - b3)

    Uz_integrand = k_tile * UR_p2 * b2 * np.exp(1j * M * theta_tile)
    Ur_integrand = (US_p2 * DrJm + (0 + 1j) * M * UT_p2 * (1.0 / r_tile) * b2) * np.exp(
        1j * M * theta_tile
    )
    Utheta_integrand = (1j * M * US_p2 * (1.0 / r_tile) * b2 - UT_p2 * DrJm) * np.exp(
        1j * M * theta_tile
    )

    Uz[:, 4] = np.real((1 / (2 * pi)) * np.trapz(Uz_integrand, k, 2))
    Ur[:, 4] = np.real((1 / (2 * pi)) * np.trapz(Ur_integrand, k, 2))
    Utheta[:, 4] = np.real((1 / (2 * pi)) * np.trapz(Utheta_integrand, k, 2))

    # okay! now they're all put together. sum up all orders and store in Greens matrices
    U_r = np.sum(Ur, 1)
    U_theta = np.sum(Utheta, 1)
    Ux = U_r * np.cos(theta) + U_theta * np.cos(theta + pi / 2)
    Uy = U_r * np.sin(theta) + U_theta * np.sin(theta + pi / 2)
    Uz = np.sum(Uz, 1)

    return Ux, Uy, Uz


def get_weights(NL, NW, L, W, nobs):
    # length
    u = np.arange(1, NL) / np.sqrt((2 * (np.arange(1, NL))) ** 2 - 1)
    vals, vecs = np.linalg.eig(np.diag(u, -1) + np.diag(u, 1))
    k = np.argsort(vals)
    bp = vals[k]
    vc = vecs[:, k]
    a, b = -L / 2, L / 2
    wL = 2 * vc[0, :] ** 2 * (b - a) / 2
    xp = (a + b) / 2 + (b - a) / 2 * bp

    # width
    u = np.arange(1, NW) / np.sqrt((2 * (np.arange(1, NW))) ** 2 - 1)
    vals, vecs = np.linalg.eig(np.diag(u, -1) + np.diag(u, 1))
    k = np.argsort(vals)
    bp = vals[k]
    vc = vecs[:, k]
    a, b = 0, W
    wW = 2 * vc[0, :] ** 2 * (b - a) / 2
    yp = (a + b) / 2 + (b - a) / 2 * bp

    # combine
    xp, yp = np.meshgrid(xp, yp)
    xp = xp.T.ravel()
    yp = yp.T.ravel()
    zp = np.zeros_like(xp)
    Xp = np.array([xp, yp, zp])

    # repmat wL and wW to the same dimensions as observation nodes for multiplication later on
    wW = np.tile(wW, (nobs, 1))
    wW = np.tile(wW[:, :, np.newaxis], (1, 1, NL))
    wL = np.tile(wL, (nobs, NW, 1))

    return wL, wW


def get_coords(NL, NW, L, W, strike, dip):
    # length
    u = np.arange(1, NL) / np.sqrt((2 * (np.arange(1, NL))) ** 2 - 1)
    vals, vecs = np.linalg.eig(np.diag(u, -1) + np.diag(u, 1))
    k = np.argsort(vals)
    bp = vals[k]
    a, b = -L / 2, L / 2
    xp = (a + b) / 2 + (b - a) / 2 * bp

    # width
    u = np.arange(1, NW) / np.sqrt((2 * (np.arange(1, NW))) ** 2 - 1)
    vals, vecs = np.linalg.eig(np.diag(u, -1) + np.diag(u, 1))
    k = np.argsort(vals)
    bp = vals[k]
    a, b = 0, W
    yp = (a + b) / 2 + (b - a) / 2 * bp

    # combine
    xp, yp = np.meshgrid(xp, yp)
    xp = xp.T.ravel()
    yp = yp.T.ravel()
    zp = np.zeros_like(xp)
    Xp = np.array([xp, yp, zp])

    R = np.array(
        [
            [1, 0, 0],
            [0, cos(np.radians(dip)), -sin(np.radians(dip))],
            [0, sin(np.radians(dip)), cos(np.radians(dip))],
        ]
    )
    Xp = R @ Xp

    R = np.array(
        [
            [cos(np.radians(strike)), -sin(np.radians(strike)), 0],
            [sin(np.radians(strike)), cos(np.radians(strike)), 0],
            [0, 0, 1],
        ]
    )
    Xp = R @ Xp
    # currently I have to reshape them BACKWARDS and then transpose them
    # I want to get .reshape(NW,NL) but because of python's indexing, that doesn't work
    xpos = Xp[0, :].reshape(NL, NW).T
    ypos = Xp[1, :].reshape(NL, NW).T
    zpos = Xp[2, :][:NW]
    return xpos, ypos, zpos

# besselj0_vectorized = np.vectorize(lambda x: mpmath.besselj(0, x))
# besselj1_vectorized = np.vectorize(lambda x: mpmath.besselj(1, x))

def backslip_cycle3D_layered_noelastic(pm, slip, xloc, T, timeEq):
    
    # reads in viscoelastic parameters from params.cfg
    consts = get_constants()

    # assumes xloc is 3, n
    # code in some kind of check for if its 2 x n or transposed
    nobs = xloc.shape[1]
    xy = xloc.copy()  

    ss, ds, ten = slip[0], slip[1], slip[2]
    kmax = 0.25
    k = np.linspace(0.000001, kmax, consts["Nterms"])
    l = pm[0]
    w = pm[1]
    dip = pm[3]
    strike = pm[4]

    dip_rad = radians(dip)
    strike_rad = radians(strike)
    d = pm[2] - w * sin(dip_rad)

    # Convert parameters to be consistent with Okada
    # if -90 <= dip <= 90:
    #    dipdir = strike_rad + pi/2
    # else:
    #    dipdir = strike_rad - pi/2

    # offset = abs(w * cos(dip_rad))
    # EastOff = offset * sin(dipdir)
    # NorthOff = offset * cos(dipdir)

    xy[0,:] -= pm[5]  # + EastOff (commented out as per MATLAB)
    xy[1,:] -= pm[6]  # + NorthOff (commented out as per MATLAB)

    Xcoord = xy[1,:]  # Switch x and y
    Ycoord = xy[0,:]

    beta1 = 1 / consts["tR1"]
    beta2 = 1 / consts["tR2"]

    # Position of point sources
    nL = int(np.ceil(1 / 40 * l))
    nW = int(np.ceil(1 / 20 * w))

    wL, wW = get_weights(nL, nW, l, w, nobs)

    # moment tensors (ss, ds, tensile)
    M1, M2, M3 = momtensor_inverse(strike_rad, dip_rad, consts["lam"], consts["Gshear"])
    # local fault coords
    xpos, ypos, zpos = get_coords(nL, nW, l, w, strike, dip)

    if ss:
        U1s = bessel_wrapper(zpos, d, consts["Gshear"], consts["lam"], M1, beta1, beta2, consts["H1"], consts["H2"], timeEq, T, k)
    if ds:
        U2s = bessel_wrapper(zpos, d, consts["Gshear"], consts["lam"], M2, beta1, beta2, consts["H1"], consts["H2"], timeEq, T, k)
    if ten:
        U3s = bessel_wrapper(zpos, d, consts["Gshear"], consts["lam"], M3, beta1, beta2, consts["H1"], consts["H2"], timeEq, T, k)

    # preallocate final greens matrices...
    G_ss = np.zeros((3, nobs, nW, nL))
    G_ds = np.zeros((3, nobs, nW, nL))
    G_ten = np.zeros((3, nobs, nW, nL))

    k_tile = np.tile(k, (nobs, 1))

    #besselj0 = np.vectorize(lambda x: float(mpmath.besselj(0, x)))
    #besselj1 = np.vectorize(lambda x: float(mpmath.besselj(1, x)))

    # we've discretized the PM up along down-dip width and along length into sub patches
    # loop through and compute displacements at surface from each of the segments
    for i, psW in enumerate(np.arange(nW)):
        for psL in np.arange(nL):

            # shift around fault
            X = Xcoord - xpos[psW, psL]
            Y = Ycoord - ypos[psW, psL]
            # compute r, theta
            r = np.sqrt(X**2 + Y**2)
            theta = np.arctan2(Y, X)

            r_tile = np.tile(r, (consts["Nterms"], 1)).T

            # calculate bessel functions
            # they're faster vectorized than using list comprehension
            # but still very slow
                    
            b0 = jv(0, (k_tile * r_tile))
            b1 = jv(1, (k_tile * r_tile))
            
            #b0 = besselj0(k_tile * r_tile)
            #b1 = besselj1(k_tile * r_tile)

            # Parallelize the computation-- takes the same time as regular scipy...
            # b0 = np.array(Parallel(n_jobs=-1)(delayed(jv)(0, x) for x in kr))
            # b1 = np.array(Parallel(n_jobs=-1)(delayed(jv)(1, x) for x in kr))

            if ss:
                G_ss[:, :, psW, psL] = integrate_disps(
                    U1s[i, :, :], k, k_tile, r_tile, theta, b0, b1
                )

            if ds:
                G_ds[:, :, psW, psL] = integrate_disps(
                    U2s[i, :, :], k, k_tile, r_tile, theta, b0, b1
                )

            if ten:
                G_ten[:, :, psW, psL] = integrate_disps(
                    U3s[i, :, :], k, k_tile, r_tile, theta, b0, b1
                )

    # for each nL x nW patch, multiply spatially based on quad weight then sum up over those patches
    # to compute gaussian quadrature integration on the pm
    G_ss = (wL * wW * G_ss).sum(axis=(2, 3))
    G_ds = (wL * wW * G_ds).sum(axis=(2, 3))
    G_ten = (wL * wW * G_ten).sum(axis=(2, 3))

    G_ds = G_ds * -1

    # solve for final disp, U
    u = ss * G_ss + ds * G_ds + ten * G_ten
    # swap the x and y back
    u = u[[1, 0, 2], :]

    u[2,:] = -1 * u[2,:]

    return u

def compute_strain(u, xs,ys):

        # Reshape U[0, :] and U[1, :] to match the grid dimensions
        Ue = u[0, :].reshape(len(xs), len(ys))
        Un = u[1, :].reshape(len(xs), len(ys))

        # Compute gradients using np.gradient
        Due_Dx, Due_Dy = np.gradient(Ue, xs, ys)  # xs and ys are the coordinates for spacing
        Dun_Dx, Dun_Dy = np.gradient(Un, xs, ys)

        # Compute strain components
        Exx = Due_Dx.ravel()
        Exy = 0.5 * (Due_Dy.ravel() + Dun_Dx.ravel())
        Eyy = Dun_Dy.ravel()

        return Exx, Exy, Eyy

def compute_strain2(u, xs,ys):

        # Reshape U[0, :] and U[1, :] to match the grid dimensions
        Ue = u[0, :].reshape(len(xs), len(ys))
        Un = u[1, :].reshape(len(xs), len(ys))

        # Compute gradients using np.gradient
        Due_Dx, Due_Dy = np.gradient(Ue, xs, ys)  # xs and ys are the coordinates for spacing
        Dun_Dx, Dun_Dy = np.gradient(Un, xs, ys)

        # Compute strain components
        Exx = Due_Dx
        Exy = 0.5 * (Due_Dy + Dun_Dx)
        Eyy = Dun_Dy

        return Exx, Exy, Eyy

# creates a regularised vector with a 20% buffer from an input vector x
# defaults to 5 km spacing
def get_vec(x, spacing=5):

    maxx = x.max()
    minx = x.min()
    dx = maxx-minx
    buffer = dx * 0.2
    n = np.rint((dx + 2 * buffer) / spacing)
    return np.linspace(minx - buffer, maxx + buffer, int(n))

def get_interseismic_strain_rates_cycles(pm, slip, x, y, T, timeEq):
    # updated interpolation method to RegularGridInterpolator as it behaves better
    # at finely meshed areas

    if np.isnan(T) | np.isnan(timeEq):

        Exx = np.zeros_like(x)
        Exy = np.zeros_like(x)
        Eyy = np.zeros_like(x)

    else:

        # from x y coordinates, build a regularized grid
        # we compute viscous deformation on a regularized grid
        # and interpolate back to our intitial input
        # xs = get_vec(x)
        # ys = get_vec(y)

        # xs = linspace(min(xystats(:,1))-100,max(xystats(:,1))+100,150);
        # ys = linspace(min(xystats(:,2)),max(xystats(:,2)),100);
        xs = np.linspace(x.min()-100,x.max()+100,150)
        ys = np.linspace(y.min(),y.max(),100)
        Xs, Ys = np.meshgrid(xs, ys)
        Xs, Ys = Xs.T.ravel(), Ys.T.ravel()
        xloc_grid = np.vstack((Xs, Ys))

        u = backslip_cycle3D_layered_noelastic(
            pm,
            slip,
            xloc_grid,
            T,
            timeEq,
        )

        Exx_ss, Exy_ss, Eyy_ss = compute_strain(u, xs, ys)
        # Exx_ss, Exy_ss, Eyy_ss = compute_strain2(u, xs, ys)
        Exx = griddata((Xs,Ys),Exx_ss,(x,y),method="linear")
        Exy = griddata((Xs,Ys),Exy_ss,(x,y),method="linear")
        Eyy = griddata((Xs,Ys),Eyy_ss,(x,y),method="linear")
        Exx = np.nan_to_num(Exx, 0)
        Exy = np.nan_to_num(Exy, 0)
        Eyy = np.nan_to_num(Eyy, 0)

        # Exx_ss, Exy_ss, Eyy_ss = compute_strain2(u, xs, ys)

        # Exx = rgt((xs,ys),Exx_ss,method='linear',bounds_error=False, fill_value = 0)((x,y))
        # Exy = rgt((xs,ys),Exy_ss,method='linear',bounds_error=False, fill_value = 0)((x,y))
        # Eyy = rgt((xs,ys),Eyy_ss,method='linear',bounds_error=False, fill_value = 0)((x,y))

    return Exx, Exy, Eyy

# add in timing...
def main():
    import numpy as np
    import matplotlib.pyplot as plt
    
    # create synthetic xy_obs for now
    xy_obs = np.random.random((2, 1000)) * 300
    xy_obs[0, :] -= 300
    xy_obs[1, :] -= 150

    # change shape
    xy_obs[0,:] = xy_obs[0,:] *2

    # hard code in a random pm, T, timeEq
    pm =np.array(
        [
            13.00682739,
            15.77467722,
            6.66666667,
            25.0,
            149.44012585,
            -363.28698341,
            116.97954325,
        ]
    )
    T = 2941.17
    timeEq = 588.235
    slip = [0,1,0]
    # Exx, Exy, Eyy = visc_tools.get_interseismic_strain_rates_cycles(pm, slip, xloc, xloc_visco_grid, T, timeEq)
    Exx, Exy, Eyy = get_interseismic_strain_rates_cycles(
        pm, slip, xy_obs[0, :], xy_obs[1, :], T, timeEq
    )

    plt.scatter(xy_obs[0,:],xy_obs[1,:],c=Eyy)
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    main()
