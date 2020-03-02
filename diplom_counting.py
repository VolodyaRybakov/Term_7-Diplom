from math import pi, sqrt, sin, cos
import scipy.constants as constants
from scipy import optimize
import time

from base_material import func_1, func_1_1, func_1_1_reversed, func_2, func_3, func_4, func_5


def findeAnomalyFromEquation(t_b: float) -> float:
    '''
        Numerical solution of the transcendental equation (by dichotomy) with respect to e_anomaly.

        returns
        -------
        e_anomaly: float'''
    e_anomaly = optimize.newton(
        func_4, pi, args=(t_b, ), maxiter=10**6)
    return e_anomaly


def checkingByAccuracy(xi_0, koef):
    time_0 = func_1_1(xi_0)
    x_s, y_s, z_s = func_1(xi_0, time_0)
    while(not func_2(x_s, y_s, z_s)):
        xi_0 += koef
        time_0 = func_1_1(xi_0)
        x_s, y_s, z_s = func_1(xi_0, time_0)
    return (xi_0 - koef, x_s, y_s, z_s)


def findTimeFromEquation(prev_t_b: float) -> float:

    xi = optimize.newton(
        func_1_1_reversed, 0, args=(prev_t_b, ), maxiter=10**6)
    return xi


if __name__ == "__main__":

    start_time = time.time()

    xi_0 = 0
    time_0 = func_1_1(xi_0)
    x_s, y_s, z_s = func_1(xi_0, time_0)
    for i in range(1, 10):
        koef = 1/(10**i)
        xi_0, x_s, y_s, z_s = checkingByAccuracy(xi_0, koef)

    time_0 = func_1_1(xi_0)
    x_s, y_s, z_s = func_1(xi_0, time_0)

    print(f"Satellite moving coords =\t({x_s},\t{y_s},\t{z_s})")

    t_0 = time_0
    print(f"t_0: {t_0}")
    t_b = func_3(t_0)

    e_anomaly_r = findeAnomalyFromEquation(t_b)
    print(f"Eccentrisity anomaly =\t\t{e_anomaly_r}")
    l_x, l_y, l_z = func_5(t_b, e_anomaly_r)
    print(f"Station angle coords =\t\t({l_x},\t{l_y},\t{l_z})")
    print(l_x**2 + l_y**2 + l_z**2)

    xi_new = findTimeFromEquation(t_b)
    print(f"Эксцентрическая аномалия:\t{xi_new}\tот времени:\t{t_b}")

    print(f"\nProgramm works {time.time() - start_time} seconds.")
