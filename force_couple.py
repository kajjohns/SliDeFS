import numpy as np
from scipy.interpolate import griddata
import warnings
from matplotlib.path import Path

def force_couple_planestrain(x1, x2, xi1, xi2, nu, mu):

    # greens function -- concentrated point force, free surface, Segall chapter 3.3
    # this solves eqns 3.36 through 3.37 for two forces across a force couple
    # first at xi1 + dx then at xi1 - dx to obtain dfdx1 for each force via first order taylor expansion
    # then do the same for dfdx2
    # xi1,xi2 are source coordinates

    # x2 is positive up
    # need to shift coordinates down, far from free surface
    x2 = x2 - 10**4
    xi2 = xi2 - 10**4

    #################################
    # dx1 couples
    dx1 = 0.001

    # first force in couple
    xi1 = xi1 + dx1

    x2minus = x2 - xi2
    x2plus = x2 + xi2
    x1minus = x1 - xi1

    r1 = np.sqrt(x1minus**2 + x2minus**2)
    r2 = np.sqrt(x1minus**2 + x2plus**2)
    theta2 = np.arctan(x1minus / x2plus)

    # line force in 1 direction
    term1 = (3 - 4 * nu) / 4 * np.log(r1)
    term2 = (8 * nu**2 - 12 * nu + 5) / 4 * np.log(r2)
    term3 = x2minus**2 / (4 * r1**2)
    term4 = ((3 - 4 * nu) * x2plus**2 + 2 * xi2 * x2plus - 2 * xi2**2) / (4 * r2**2)
    term5 = -xi2 * x2 * x2plus**2 / r2**4
    g11_1 = (
        -1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4 + term5)
    )  # eqn. 3.36

    term1 = (1 - 2 * nu) * (1 - nu) * theta2
    term2 = x2minus * x1minus / (4 * r1**2)
    term3 = (3 - 4 * nu) * x2minus * x1minus / (4 * r2**2)
    term4 = -xi2 * x2 * x1minus * x2plus / r2**4
    g21_1 = (
        1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4)
    )  # eqn 3.37

    term1 = x1minus**3 / r1**4
    term2 = x1minus * (x1minus**2 - 4 * xi2 * x2 - 2 * xi2**2) / r2**4
    term3 = 8 * xi2 * x2 * x1minus * x2plus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x1minus / r1**2 + 3 * x1minus / r2**2 - 4 * x2 * x1minus * x2plus / r2**4)
    )
    s11_1_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)  # eqn 3.42

    term1 = x2minus * x1minus**2 / r1**4
    term2 = x2plus * (2 * xi2 * x2 + x1minus**2) / r2**4
    term3 = -8 * xi2 * x2 * x1minus**2 * x2plus / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x2minus / r1**2 + (3 * x2 + xi2) / r2**2 - 4 * x2 * x2plus**2 / r2**4)
    )
    s12_1_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)  # eqn 3.43

    term1 = x2minus**2 * x1minus / r1**4
    term2 = -x1minus * (xi2**2 - x2**2 + 6 * xi2 * x2) / r2**4
    term3 = 8 * xi2 * x2 * x1minus**3 / r2**6
    term4 = (
        -(1 - 2 * nu)
        / 2
        * (x1minus / r1**2 - x1minus / r2**2 - 4 * x2 * x1minus * x2plus / r2**4)
    )
    s22_1_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)  # eqn 3.44

    # line force in 2 direction
    term1 = -(1 - 2 * nu) * (1 - nu) * theta2
    term2 = x2minus * x1minus / (4 * r1**2)
    term3 = (3 - 4 * nu) * x2minus * x1minus / (4 * r2**2)
    term4 = xi2 * x2 * x1minus * x2plus / r2**4
    g12_1 = (
        1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4)
    )  # eqn 3.38

    term1 = -(3 - 4 * nu) / 4 * np.log(r1)
    term2 = -(8 * nu**2 - 12 * nu + 5) / 4 * np.log(r2)
    term3 = -(x1minus**2) / (4 * r1**2)
    term4 = (2 * xi2 * x2 - (3 - 4 * nu) * x1minus**2) / (4 * r2**2)
    term5 = -xi2 * x2 * x1minus**2 / r2**4
    g22_1 = (
        1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4 + term5)
    )  # eqn 3.39

    term1 = x2minus * x1minus**2 / r1**4
    term2 = (x2plus * (x1minus**2 + 2 * xi2**2) - 2 * xi2 * x1minus**2) / r2**4
    term3 = 8 * xi2 * x2 * x2plus * x1minus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (-x2minus / r1**2 + (3 * xi2 + x2) / r2**2 + 4 * x2 * x1minus**2 / r2**4)
    )
    s11_2_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)  # eqn 3.45

    term1 = x2minus**2 / r1**4
    term2 = (x2**2 - 2 * xi2 * x2 - xi2**2) / r2**4
    term3 = 8 * xi2 * x2 * x2plus**2 / r2**6
    term4 = (1 - 2 * nu) / 2 * (1 / r1**2 - 1 / r2**2 + 4 * x2 * x2plus / r2**4)
    s12_2_1 = (
        -x1minus / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    )  # eqn 3.46

    term1 = x2minus**3 / r1**4
    term2 = x2plus * (x2plus**2 + 2 * xi2 * x2) / r2**4
    term3 = -8 * xi2 * x2 * x2plus * x1minus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x2minus / r1**2 + (3 * x2 + xi2) / r2**2 - 4 * x2 * x1minus**2 / r2**4)
    )
    s22_2_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)  # eqn 3.47

    # second force in couple

    xi1 = xi1 - dx1 - dx1

    x2minus = x2 - xi2
    x2plus = x2 + xi2
    x1minus = x1 - xi1

    r1 = np.sqrt(x1minus**2 + x2minus**2)
    r2 = np.sqrt(x1minus**2 + x2plus**2)
    theta2 = np.arctan(x1minus / x2plus)

    # line force in 1 direction
    term1 = (3 - 4 * nu) / 4 * np.log(r1)
    term2 = (8 * nu**2 - 12 * nu + 5) / 4 * np.log(r2)
    term3 = x2minus**2 / (4 * r1**2)
    term4 = ((3 - 4 * nu) * x2plus**2 + 2 * xi2 * x2plus - 2 * xi2**2) / (4 * r2**2)
    term5 = -xi2 * x2 * x2plus**2 / r2**4
    g11_2 = -1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4 + term5)
    term1 = (1 - 2 * nu) * (1 - nu) * theta2
    term2 = x2minus * x1minus / (4 * r1**2)
    term3 = (3 - 4 * nu) * x2minus * x1minus / (4 * r2**2)
    term4 = -xi2 * x2 * x1minus * x2plus / r2**4
    g21_2 = 1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4)

    term1 = x1minus**3 / r1**4
    term2 = x1minus * (x1minus**2 - 4 * xi2 * x2 - 2 * xi2**2) / r2**4
    term3 = 8 * xi2 * x2 * x1minus * x2plus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x1minus / r1**2 + 3 * x1minus / r2**2 - 4 * x2 * x1minus * x2plus / r2**4)
    )
    s11_1_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus * x1minus**2 / r1**4
    term2 = x2plus * (2 * xi2 * x2 + x1minus**2) / r2**4
    term3 = -8 * xi2 * x2 * x1minus**2 * x2plus / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x2minus / r1**2 + (3 * x2 + xi2) / r2**2 - 4 * x2 * x2plus**2 / r2**4)
    )
    s12_1_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus**2 * x1minus / r1**4
    term2 = -x1minus * (xi2**2 - x2**2 + 6 * xi2 * x2) / r2**4
    term3 = 8 * xi2 * x2 * x1minus**3 / r2**6
    term4 = (
        -(1 - 2 * nu)
        / 2
        * (x1minus / r1**2 - x1minus / r2**2 - 4 * x2 * x1minus * x2plus / r2**4)
    )
    s22_1_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)

    # line force in 2 direction
    term1 = -(1 - 2 * nu) * (1 - nu) * theta2
    term2 = x2minus * x1minus / (4 * r1**2)
    term3 = (3 - 4 * nu) * x2minus * x1minus / (4 * r2**2)
    term4 = xi2 * x2 * x1minus * x2plus / r2**4
    g12_2 = 1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = -(3 - 4 * nu) / 4 * np.log(r1)
    term2 = -(8 * nu**2 - 12 * nu + 5) / 4 * np.log(r2)
    term3 = -(x1minus**2) / (4 * r1**2)
    term4 = (2 * xi2 * x2 - (3 - 4 * nu) * x1minus**2) / (4 * r2**2)
    term5 = -xi2 * x2 * x1minus**2 / r2**4
    g22_2 = 1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4 + term5)

    term1 = x2minus * x1minus**2 / r1**4
    term2 = (x2plus * (x1minus**2 + 2 * xi2**2) - 2 * xi2 * x1minus**2) / r2**4
    term3 = 8 * xi2 * x2 * x2plus * x1minus**2 / r2**6
    term4 = ((1 - 2 * nu)  / 2  * (-x2minus / r1**2 + (3 * xi2 + x2) / r2**2 + 4 * x2 * x1minus**2 / r2**4))
    s11_2_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus**2 / r1**4
    term2 = (x2**2 - 2 * xi2 * x2 - xi2**2) / r2**4
    term3 = 8 * xi2 * x2 * x2plus**2 / r2**6
    term4 = (1 - 2 * nu) / 2 * (1 / r1**2 - 1 / r2**2 + 4 * x2 * x2plus / r2**4)
    s12_2_2 = -x1minus / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus**3 / r1**4
    term2 = x2plus * (x2plus**2 + 2 * xi2 * x2) / r2**4
    term3 = -8 * xi2 * x2 * x2plus * x1minus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x2minus / r1**2 + (3 * x2 + xi2) / r2**2 - 4 * x2 * x1minus**2 / r2**4)
    )
    s22_2_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)

    df1dx = {}
    df1dx["s11"] = (s11_1_1 - s11_1_2) / (2 * dx1)
    df1dx["s12"] = (s12_1_1 - s12_1_2) / (2 * dx1)
    df1dx["s22"] = (s22_1_1 - s22_1_2) / (2 * dx1)
    df1dx["u1"] = (g11_1 - g11_2) / (2 * dx1)
    df1dx["u2"] = (g21_1 - g21_2) / (2 * dx1)

    df2dx = {}
    df2dx["s11"] = (s11_2_1 - s11_2_2) / (2 * dx1)
    df2dx["s12"] = (s12_2_1 - s12_2_2) / (2 * dx1)
    df2dx["s22"] = (s22_2_1 - s22_2_2) / (2 * dx1)
    df2dx["u1"] = (g12_1 - g12_2) / (2 * dx1)
    df2dx["u2"] = (g22_1 - g22_2) / (2 * dx1)

    xi1 = xi1 + dx1 + dx1

    #################################
    # dx2 couples
    dx2 = 0.001

    # first force in couple
    xi2 = xi2 + dx2

    x2minus = x2 - xi2
    x2plus = x2 + xi2
    x1minus = x1 - xi1

    r1 = np.sqrt(x1minus**2 + x2minus**2)
    r2 = np.sqrt(x1minus**2 + x2plus**2)
    theta2 = np.arctan(x1minus / x2plus)

    # line force in 1 direction
    term1 = (3 - 4 * nu) / 4 * np.log(r1)
    term2 = (8 * nu**2 - 12 * nu + 5) / 4 * np.log(r2)
    term3 = x2minus**2 / (4 * r1**2)
    term4 = ((3 - 4 * nu) * x2plus**2 + 2 * xi2 * x2plus - 2 * xi2**2) / (4 * r2**2)
    term5 = -xi2 * x2 * x2plus**2 / r2**4
    g11_1 = -1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4 + term5)
    term1 = (1 - 2 * nu) * (1 - nu) * theta2
    term2 = x2minus * x1minus / (4 * r1**2)
    term3 = (3 - 4 * nu) * x2minus * x1minus / (4 * r2**2)
    term4 = -xi2 * x2 * x1minus * x2plus / r2**4
    g21_1 = 1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4)

    term1 = x1minus**3 / r1**4
    term2 = x1minus * (x1minus**2 - 4 * xi2 * x2 - 2 * xi2**2) / r2**4
    term3 = 8 * xi2 * x2 * x1minus * x2plus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x1minus / r1**2 + 3 * x1minus / r2**2 - 4 * x2 * x1minus * x2plus / r2**4)
    )
    s11_1_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus * x1minus**2 / r1**4
    term2 = x2plus * (2 * xi2 * x2 + x1minus**2) / r2**4
    term3 = -8 * xi2 * x2 * x1minus**2 * x2plus / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x2minus / r1**2 + (3 * x2 + xi2) / r2**2 - 4 * x2 * x2plus**2 / r2**4)
    )
    s12_1_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus**2 * x1minus / r1**4
    term2 = -x1minus * (xi2**2 - x2**2 + 6 * xi2 * x2) / r2**4
    term3 = 8 * xi2 * x2 * x1minus**3 / r2**6
    term4 = (
        -(1 - 2 * nu)
        / 2
        * (x1minus / r1**2 - x1minus / r2**2 - 4 * x2 * x1minus * x2plus / r2**4)
    )
    s22_1_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)

    # line force in 2 direction
    term1 = -(1 - 2 * nu) * (1 - nu) * theta2
    term2 = x2minus * x1minus / (4 * r1**2)
    term3 = (3 - 4 * nu) * x2minus * x1minus / (4 * r2**2)
    term4 = xi2 * x2 * x1minus * x2plus / r2**4
    g12_1 = 1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = -(3 - 4 * nu) / 4 * np.log(r1)
    term2 = -(8 * nu**2 - 12 * nu + 5) / 4 * np.log(r2)
    term3 = -(x1minus**2) / (4 * r1**2)
    term4 = (2 * xi2 * x2 - (3 - 4 * nu) * x1minus**2) / (4 * r2**2)
    term5 = -xi2 * x2 * x1minus**2 / r2**4
    g22_1 = 1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4 + term5)

    term1 = x2minus * x1minus**2 / r1**4
    term2 = (x2plus * (x1minus**2 + 2 * xi2**2) - 2 * xi2 * x1minus**2) / r2**4
    term3 = 8 * xi2 * x2 * x2plus * x1minus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (-x2minus / r1**2 + (3 * xi2 + x2) / r2**2 + 4 * x2 * x1minus**2 / r2**4)
    )
    s11_2_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus**2 / r1**4
    term2 = (x2**2 - 2 * xi2 * x2 - xi2**2) / r2**4
    term3 = 8 * xi2 * x2 * x2plus**2 / r2**6
    term4 = (1 - 2 * nu) / 2 * (1 / r1**2 - 1 / r2**2 + 4 * x2 * x2plus / r2**4)
    s12_2_1 = -x1minus / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus**3 / r1**4
    term2 = x2plus * (x2plus**2 + 2 * xi2 * x2) / r2**4
    term3 = -8 * xi2 * x2 * x2plus * x1minus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x2minus / r1**2 + (3 * x2 + xi2) / r2**2 - 4 * x2 * x1minus**2 / r2**4)
    )
    s22_2_1 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)

    # second force in couple

    xi2 = xi2 - dx2 - dx2

    x2minus = x2 - xi2
    x2plus = x2 + xi2
    x1minus = x1 - xi1

    r1 = np.sqrt(x1minus**2 + x2minus**2)
    r2 = np.sqrt(x1minus**2 + x2plus**2)
    theta2 = np.arctan(x1minus / x2plus)

    # line force in 1 direction
    term1 = (3 - 4 * nu) / 4 * np.log(r1)
    term2 = (8 * nu**2 - 12 * nu + 5) / 4 * np.log(r2)
    term3 = x2minus**2 / (4 * r1**2)
    term4 = ((3 - 4 * nu) * x2plus**2 + 2 * xi2 * x2plus - 2 * xi2**2) / (4 * r2**2)
    term5 = -xi2 * x2 * x2plus**2 / r2**4
    g11_2 = -1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4 + term5)
    term1 = (1 - 2 * nu) * (1 - nu) * theta2
    term2 = x2minus * x1minus / (4 * r1**2)
    term3 = (3 - 4 * nu) * x2minus * x1minus / (4 * r2**2)
    term4 = -xi2 * x2 * x1minus * x2plus / r2**4
    g21_2 = 1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4)

    term1 = x1minus**3 / r1**4
    term2 = x1minus * (x1minus**2 - 4 * xi2 * x2 - 2 * xi2**2) / r2**4
    term3 = 8 * xi2 * x2 * x1minus * x2plus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x1minus / r1**2 + 3 * x1minus / r2**2 - 4 * x2 * x1minus * x2plus / r2**4)
    )
    s11_1_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus * x1minus**2 / r1**4
    term2 = x2plus * (2 * xi2 * x2 + x1minus**2) / r2**4
    term3 = -8 * xi2 * x2 * x1minus**2 * x2plus / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x2minus / r1**2 + (3 * x2 + xi2) / r2**2 - 4 * x2 * x2plus**2 / r2**4)
    )
    s12_1_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus**2 * x1minus / r1**4
    term2 = -x1minus * (xi2**2 - x2**2 + 6 * xi2 * x2) / r2**4
    term3 = 8 * xi2 * x2 * x1minus**3 / r2**6
    term4 = (
        -(1 - 2 * nu)
        / 2
        * (x1minus / r1**2 - x1minus / r2**2 - 4 * x2 * x1minus * x2plus / r2**4)
    )
    s22_1_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)

    # line force in 2 direction
    term1 = -(1 - 2 * nu) * (1 - nu) * theta2
    term2 = x2minus * x1minus / (4 * r1**2)
    term3 = (3 - 4 * nu) * x2minus * x1minus / (4 * r2**2)
    term4 = xi2 * x2 * x1minus * x2plus / r2**4
    g12_2 = 1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = -(3 - 4 * nu) / 4 * np.log(r1)
    term2 = -(8 * nu**2 - 12 * nu + 5) / 4 * np.log(r2)
    term3 = -(x1minus**2) / (4 * r1**2)
    term4 = (2 * xi2 * x2 - (3 - 4 * nu) * x1minus**2) / (4 * r2**2)
    term5 = -xi2 * x2 * x1minus**2 / r2**4
    g22_2 = 1 / (2 * np.pi * mu * (1 - nu)) * (term1 + term2 + term3 + term4 + term5)

    term1 = x2minus * x1minus**2 / r1**4
    term2 = (x2plus * (x1minus**2 + 2 * xi2**2) - 2 * xi2 * x1minus**2) / r2**4
    term3 = 8 * xi2 * x2 * x2plus * x1minus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (-x2minus / r1**2 + (3 * xi2 + x2) / r2**2 + 4 * x2 * x1minus**2 / r2**4)
    )
    s11_2_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus**2 / r1**4
    term2 = (x2**2 - 2 * xi2 * x2 - xi2**2) / r2**4
    term3 = 8 * xi2 * x2 * x2plus**2 / r2**6
    term4 = (1 - 2 * nu) / 2 * (1 / r1**2 - 1 / r2**2 + 4 * x2 * x2plus / r2**4)
    s12_2_2 = -x1minus / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)
    term1 = x2minus**3 / r1**4
    term2 = x2plus * (x2plus**2 + 2 * xi2 * x2) / r2**4
    term3 = -8 * xi2 * x2 * x2plus * x1minus**2 / r2**6
    term4 = (
        (1 - 2 * nu)
        / 2
        * (x2minus / r1**2 + (3 * x2 + xi2) / r2**2 - 4 * x2 * x1minus**2 / r2**4)
    )
    s22_2_2 = -1 / (2 * np.pi * (1 - nu)) * (term1 + term2 + term3 + term4)

    df1dz = {}
    df1dz["s11"] = (s11_1_1 - s11_1_2) / (2 * dx2)
    df1dz["s12"] = (s12_1_1 - s12_1_2) / (2 * dx2)
    df1dz["s22"] = (s22_1_1 - s22_1_2) / (2 * dx2)
    df1dz["u1"] = (g11_1 - g11_2) / (2 * dx2)
    df1dz["u2"] = (g21_1 - g21_2) / (2 * dx2)

    df2dz = {}
    df2dz["s11"] = (s11_2_1 - s11_2_2) / (2 * dx2)
    df2dz["s12"] = (s12_2_1 - s12_2_2) / (2 * dx2)
    df2dz["s22"] = (s22_2_1 - s22_2_2) / (2 * dx2)
    df2dz["u1"] = (g12_1 - g12_2) / (2 * dx2)
    df2dz["u2"] = (g22_1 - g22_2) / (2 * dx2)

    return df1dx, df2dx, df1dz, df2dz


# This next series of functions is used in the triangular moment source calculations
def make_grid(S, T, xobs, yobs):

    ns = len(S)
    nt = len(T)
    nobs = len(xobs)

    Sn = np.tile(S, (nt, 1))
    Tn = np.tile(T, (ns, 1)).T

    Sn = np.tile(Sn, (nobs, 1, 1))
    Tn = np.tile(Tn, (nobs, 1, 1))

    xn = np.tile(xobs, (nt, ns, 1)).transpose(2, 0, 1)
    yn = np.tile(yobs, (nt, ns, 1)).transpose(2, 0, 1)

    return Sn, Tn, xn, yn

def calc_dg11dx1(xo, yo, nd, s, t, nu):
    
    t1 = (
        2
        * (1 + nu)
        * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
        * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
    )
    t2 = (
        (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
    ) ** 2
    t3 = (3 - nu) * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
    t4 = (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2 + (
        -((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo
    ) ** 2
    return t1 / t2 - t3 / t4

def calc_dg11dx2(xo, yo, nd, s, t, nu):

    return (
        (
            2
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 3
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        - (
            (3 - nu)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        - (
            2
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
    )

def calc_dg12dx1(xo, yo, nd, s, t, nu):
    return (
        2
        * (-1 - nu)
        * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
    ) / (
        (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
    ) ** 2 - (
        (-1 - nu) * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
    ) / (
        (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
    )

def calc_dg12dx2(xo, yo, nd, s, t, nu):
    return (
        2
        * (-1 - nu)
        * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
        * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
    ) / (
        (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
    ) ** 2 - (
        (-1 - nu) * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
    ) / (
        (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
    )

def calc_dg22dx1(xo, yo, nd, s, t, nu):
    return (
        (
            2
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 3
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        - (
            (3 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        - (
            2
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
    )

def calc_dg22dx2(xo, yo, nd, s, t, nu):
    return (
        2
        * (1 + nu)
        * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
    ) / (
        (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
    ) ** 2 - (
        (3 - nu) * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
    ) / (
        (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
    )

def calc_Ddg11dx1Dx(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (1 + nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                ** 2
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (3 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            2
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        - (3 - nu)
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
    )

def calc_Ddg11dx1Dy(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (1 + nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 3
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (3 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            4
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
    )

def calc_Ddg11dx2Dx(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (1 + nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 3
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (3 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            4
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
    )

def calc_Ddg11dx2Dy(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (1 + nu)
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 4
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (3 - nu)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            10
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        - (3 - nu)
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        - (2 * (1 + nu))
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
    )

def calc_Ddg12dx1Dx(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (-1 - nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                ** 3
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            6
            * (-1 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
    )

def calc_Ddg12dx1Dy(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (-1 - nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                ** 2
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (-1 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            2
            * (-1 - nu)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        - (-1 - nu)
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
    )

def calc_Ddg12dx2Dx(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (-1 - nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                ** 2
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (-1 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            2
            * (-1 - nu)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        - (-1 - nu)
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
    )

def calc_Ddg12dx2Dy(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (-1 - nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 3
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            6
            * (-1 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
    )

def calc_Ddg22dx1Dx(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (1 + nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                ** 4
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (3 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            10
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        - (3 - nu)
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        - (2 * (1 + nu))
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
    )

def calc_Ddg22dx1Dy(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (1 + nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                ** 3
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (3 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            4
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
    )

def calc_Ddg22dx2Dx(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (1 + nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                ** 3
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (3 - nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            4
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
    )

def calc_Ddg22dx2Dy(xo, yo, nd, s, t, nu):
    return (
        -(
            (
                8
                * (1 + nu)
                * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo)
                ** 2
                * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            / (
                (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
                + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo)
                ** 2
            )
            ** 3
        )
        + (
            2
            * (1 + nu)
            * (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        + (
            2
            * (3 - nu)
            * (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
        ** 2
        - (3 - nu)
        / (
            (-((1 - s - t) * nd[0, 0]) - t * nd[1, 0] - s * nd[2, 0] + xo) ** 2
            + (-((1 - s - t) * nd[0, 1]) - t * nd[1, 1] - s * nd[2, 1] + yo) ** 2
        )
    )


def triangular_moment_source(nd, obs_xyz, nu):

    xo, yo = obs_xyz[:, :2].T

    dx = (nd[:, 0].max() - nd[:, 0].min()) / 10
    dy = (nd[:, 1].max() - nd[:, 1].min()) / 10

    centroid_xy = np.array(
        [
            [nd[:, 0].mean() + dx / 2, nd[:, 1].mean()],
            [nd[:, 0].mean() - dx / 2, nd[:, 1].mean()],
            [nd[:, 0].mean(), nd[:, 1].mean() + dy / 2],
            [nd[:, 0].mean(), nd[:, 1].mean() - dy / 2],
        ]
    )

    # Integration grid
    T = np.linspace(0, 1, 10)
    S = np.linspace(0, 1, 10)
    ds, dt = S[1] - S[0], T[1] - T[0]

    T_grid, S_grid = np.meshgrid(T, S, indexing="ij")
    mask = S_grid <= (1 - T_grid)
    # Create grids for centroid-based calculations
    Sv, Tv, xv, yv = make_grid(S, T, centroid_xy[:, 0], centroid_xy[:, 1])

    # build list of components as STR which correspond to calc_STR() functions below
    # for legibility, we loop through each component and assign the variable STR =  calc_STR(args)
    #
    disp_comps = [
        f"dg{comp}d{dir}" for dir in ["x1", "x2"] for comp in ["11", "12", "22"]
    ]
    strain_comps = [f"D{fn}{d}" for fn in disp_comps for d in ["Dx", "Dy"]]
    all_comps = disp_comps + strain_comps
    results = {}
    # Compute displacement components at centroids
    for fn in disp_comps:
        results[f"{fn}_centroid"] = -(ds * dt) * globals()[f"calc_{fn}"](
            xv, yv, nd, Sv, Tv, nu
        )[:, mask].sum(axis=1)

    # Create grids for displacement and strain calculations
    Sv, Tv, xv, yv = make_grid(S, T, xo, yo)

    # Compute displacement and strain components
    for fn in all_comps:
        results[fn] = -(ds * dt) * globals()[f"calc_{fn}"](xv, yv, nd, Sv, Tv, nu)[
            :, mask
        ].sum(axis=1)

    # for points inside the node, replace computed strain componentes
    # with finite difference computed using centroid displacements
    polygon_path = Path(nd)
    inside = polygon_path.contains_points(obs_xyz[:, :2])

    indices = [(0, 1, "Dx", dx), (2, 3, "Dy", dy)]
    for fn in disp_comps:
        for i1, i2, d_suffix, d_var in indices:
            results[f"D{fn}{d_suffix}"][inside] = (
                results[f"{fn}_centroid"][i1] - results[f"{fn}_centroid"][i2]
            ) / d_var
    # Assign symmetrical components
    for key in [
        "dg21dx1",
        "dg21dx2",
        "Ddg21dx1Dx",
        "Ddg21dx1Dy",
        "Ddg21dx2Dx",
        "Ddg21dx2Dy",
    ]:
        results[key] = results[key.replace("21", "12")]

    def build_matrix(keys):
        return {key: results[val] for key, val in keys.items()}

    u1 = build_matrix(
        {"m11": "dg11dx1", "m12": "dg11dx2", "m21": "dg12dx1", "m22": "dg12dx2"}
    )

    u2 = build_matrix(
        {"m11": "dg21dx1", "m12": "dg21dx2", "m21": "dg22dx1", "m22": "dg22dx2"}
    )

    e11 = build_matrix(
        {
            "m11": "Ddg11dx1Dx",
            "m12": "Ddg11dx2Dx",
            "m21": "Ddg12dx1Dx",
            "m22": "Ddg12dx2Dx",
        }
    )

    e22 = build_matrix(
        {
            "m11": "Ddg21dx1Dy",
            "m12": "Ddg21dx2Dy",
            "m21": "Ddg22dx1Dy",
            "m22": "Ddg22dx2Dy",
        }
    )

    e12 = {
        "m11": 0.5 * (results["Ddg11dx1Dy"] + results["Ddg21dx1Dx"]),
        "m12": 0.5 * (results["Ddg11dx2Dy"] + results["Ddg21dx2Dx"]),
        "m21": 0.5 * (results["Ddg12dx1Dy"] + results["Ddg22dx1Dx"]),
        "m22": 0.5 * (results["Ddg12dx2Dy"] + results["Ddg22dx2Dx"]),
    }

    return u1, u2, e11, e12, e22


def buildG_MomentSource_2d(mesh, xy_obs, nu_0):
    # returns moment sources computed as points

    GExx_m11 = np.zeros((len(xy_obs), len(mesh["tri_centroids"])))
    GExx_m12 = np.zeros_like(GExx_m11)
    GExx_m22 = np.zeros_like(GExx_m11)

    GExy_m11 = np.zeros_like(GExx_m11)
    GExy_m12 = np.zeros_like(GExx_m11)
    GExy_m22 = np.zeros_like(GExx_m11)

    GEyy_m11 = np.zeros_like(GExx_m11)
    GEyy_m12 = np.zeros_like(GExx_m11)
    GEyy_m22 = np.zeros_like(GExx_m11)

    # convert from plane strain to plane stress by transforming elastic
    # properties
    # builds moment sources as numerical integrations over triangular elements
    nu = nu_0 / (
        1 + nu_0
    )  # Note: no need to transform Lame constants because solution only depends on nu
    # GFs depend only on ratios of Lame constants
    mu = 1.0
    lam = 2 * nu / (1 - 2 * nu)

    for k, (x_i, y_i) in enumerate(mesh["nodes"][:, 0:2]):
        # compute response at centroids due to a point source at nodes to avoid singularities
        [df1dx, df2dx, df1dy, df2dy] = force_couple_planestrain(
            mesh["tri_centroids"][:, 0], mesh["tri_centroids"][:, 1], x_i, y_i, nu, mu
        )

        # force_couple_planesrain returns stresses - need to to convert to strain
        # rate.

        df1dx["exx"] = ((lam + 2 * mu) * df1dx["s11"] - lam * df1dx["s22"]) / (
            2 * mu * (lam + 2 * mu) + lam * (lam + 2 * mu) - lam**2
        )
        df1dx["exy"] = df1dx["s12"] / (2 * mu)
        df1dx["eyy"] = (df1dx["s22"] - lam * df1dx["exx"]) / (lam + 2 * mu)

        df2dx["exx"] = ((lam + 2 * mu) * df2dx["s11"] - lam * df2dx["s22"]) / (
            2 * mu * (lam + 2 * mu) + lam * (lam + 2 * mu) - lam**2
        )
        df2dx["exy"] = df2dx["s12"] / (2 * mu)
        df2dx["eyy"] = (df2dx["s22"] - lam * df2dx["exx"]) / (lam + 2 * mu)

        df1dy["exx"] = ((lam + 2 * mu) * df1dy["s11"] - lam * df1dy["s22"]) / (
            2 * mu * (lam + 2 * mu) + lam * (lam + 2 * mu) - lam**2
        )
        df1dy["exy"] = df1dy["s12"] / (2 * mu)
        df1dy["eyy"] = (df1dy["s22"] - lam * df1dy["exx"]) / (lam + 2 * mu)

        df2dy["exx"] = ((lam + 2 * mu) * df2dy["s11"] - lam * df2dy["s22"]) / (
            2 * mu * (lam + 2 * mu) + lam * (lam + 2 * mu) - lam**2
        )
        df2dy["exy"] = df2dy["s12"] / (2 * mu)
        df2dy["eyy"] = (df2dy["s22"] - lam * df2dy["exx"]) / (lam + 2 * mu)

        exx_m11 = df1dx["exx"]
        exx_m12 = df1dy["exx"] + df2dx["exx"]
        exx_m22 = df2dy["exx"]

        exy_m11 = df1dx["exy"]
        exy_m12 = df1dy["exy"] + df2dx["exy"]
        exy_m22 = df2dy["exy"]

        eyy_m11 = df1dx["eyy"]
        eyy_m12 = df1dy["eyy"] + df2dx["eyy"]
        eyy_m22 = df2dy["eyy"]

        # reinterpolate to desired grid
        GExx_m11[:, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exx_m11, xy_obs[:, 0:2], method="linear"
        )
        GExx_m12[:, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exx_m12, xy_obs[:, 0:2], method="linear"
        )
        GExx_m22[:, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exx_m22, xy_obs[:, 0:2], method="linear"
        )

        idx = np.isnan(GExx_m11[:, k])
        GExx_m11[idx, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exx_m11, xy_obs[idx, 0:2], method="nearest"
        )
        GExx_m12[idx, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exx_m12, xy_obs[idx, 0:2], method="nearest"
        )
        GExx_m22[idx, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exx_m22, xy_obs[idx, 0:2], method="nearest"
        )

        GExy_m11[:, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exy_m11, xy_obs[:, 0:2], method="linear"
        )
        GExy_m12[:, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exy_m12, xy_obs[:, 0:2], method="linear"
        )
        GExy_m22[:, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exy_m22, xy_obs[:, 0:2], method="linear"
        )

        GExy_m11[idx, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exy_m11, xy_obs[idx, 0:2], method="nearest"
        )
        GExy_m12[idx, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exy_m12, xy_obs[idx, 0:2], method="nearest"
        )
        GExy_m22[idx, k] = griddata(
            mesh["tri_centroids"][:, 0:2], exy_m22, xy_obs[idx, 0:2], method="nearest"
        )

        GEyy_m11[:, k] = griddata(
            mesh["tri_centroids"][:, 0:2], eyy_m11, xy_obs[:, 0:2], method="linear"
        )
        GEyy_m12[:, k] = griddata(
            mesh["tri_centroids"][:, 0:2], eyy_m12, xy_obs[:, 0:2], method="linear"
        )
        GEyy_m22[:, k] = griddata(
            mesh["tri_centroids"][:, 0:2], eyy_m22, xy_obs[:, 0:2], method="linear"
        )

        GEyy_m11[idx, k] = griddata(
            mesh["tri_centroids"][:, 0:2], eyy_m11, xy_obs[idx, 0:2], method="nearest"
        )
        GEyy_m12[idx, k] = griddata(
            mesh["tri_centroids"][:, 0:2], eyy_m12, xy_obs[idx, 0:2], method="nearest"
        )
        GEyy_m22[idx, k] = griddata(
            mesh["tri_centroids"][:, 0:2], eyy_m22, xy_obs[idx, 0:2], method="nearest"
        )

    GExx = {"m11": GExx_m11, "m12": GExx_m12, "m22": GExx_m22}
    GExy = {"m11": GExy_m11, "m12": GExy_m12, "m22": GExy_m22}
    GEyy = {"m11": GEyy_m11, "m12": GEyy_m12, "m22": GEyy_m22}
    G = {"Exx": GExx, "Exy": GExy, "Eyy": GEyy}
    return G


def buildG_MomentSource_2d_tri(mesh, obs_xyz, nu_0):
    # builds moment sources as numerical integrations over triangular elements
    nu = nu_0 / (1 + nu_0)
    # Preallocate arrays
    G_mom = {
        f"G{comp}_{mom}": np.zeros((len(obs_xyz), len(mesh["tri"])))
        for comp in ["Exx", "Exy", "Eyy"]
        for mom in ["m11", "m12", "m22"]
    }

    warnings.filterwarnings("ignore")
    for k, nd in enumerate(mesh["elts"]):
        # first node to feed into triangular moment source
        u1, u2, e11, e12, e22 = triangular_moment_source(
            nd, mesh["tri_centroids"][:, :2], nu
        )
        e11["m12"] = e11["m12"] + e11["m21"]
        e12["m12"] = e12["m12"] + e12["m21"]
        e22["m12"] = e22["m12"] + e22["m21"]

        Es = {"Exx": e11, "Exy": e12, "Eyy": e22}
        moms = ["m11", "m12", "m22"]

        for comp in Es:
            for mom in moms:
                Gkey = f"G{comp}_{mom}"
                data = Es[comp][mom]
                G_mom[Gkey][:, k] = griddata(
                    mesh["tri_centroids"][:, :2], data, obs_xyz[:, :2], method="linear"
                )
                idx = np.isnan(G_mom[Gkey][:, k])
                G_mom[Gkey][idx, k] = griddata(
                    mesh["tri_centroids"][:, :2], data, obs_xyz[idx, :2], method="nearest"
                )


    warnings.filterwarnings("default")

    GExx = {
        "m11": G_mom["GExx_m11"],
        "m12": G_mom["GExx_m12"],
        "m22": G_mom["GExx_m22"],
    }
    GExy = {
        "m11": G_mom["GExy_m11"],
        "m12": G_mom["GExy_m12"],
        "m22": G_mom["GExy_m22"],
    }
    GEyy = {
        "m11": G_mom["GEyy_m11"],
        "m12": G_mom["GEyy_m12"],
        "m22": G_mom["GEyy_m22"],
    }

    G = {"Exx": GExx, "Exy": GExy, "Eyy": GEyy}
    return G
