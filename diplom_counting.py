from math import pi, sqrt, sin, cos
import scipy.constants as constants
from scipy import optimize
import time

from base_material import func_1, func_1_1, func_1_1_reversed, func_2, func_3, func_4, func_5


def checkingByAccuracy(xi_0, koef):
    time_0 = func_1_1(xi_0)
    x_s, y_s, z_s = func_1(xi_0, time_0)
    while(not func_2(x_s, y_s, z_s)):
        xi_0 += koef
        time_0 = func_1_1(xi_0)
        x_s, y_s, z_s = func_1(xi_0, time_0)
    return xi_0 - koef


def getEAnomalyForFirstLocation(xi_0):
    for i in range(1, 10):
        koef = 1/(10**i)
        xi_0 = checkingByAccuracy(xi_0, koef)

    return xi_0


def getFirstLocationAngleCoords():
    '''
        Задача 1
        Найти начальное время и начальный угол наведения телескопа'''

    xi_0 = 0

    xi_0 = getEAnomalyForFirstLocation(xi_0)

    time_0 = func_1_1(xi_0)
    x_s, y_s, z_s = func_1(xi_0, time_0)

    print(f"Satellite moving coords =\t({x_s},\t{y_s},\t{z_s})")

    t_0 = time_0
    print(f"t_0: {t_0}")
    t_b = func_3(t_0)

    xi_r = func_4(t_b)
    print(f"Eccentrisity anomaly =\t\t{xi_r}")
    l_x, l_y, l_z = func_5(t_b, xi_r)
    print(f"Station angle coords =\t\t({l_x},\t{l_y},\t{l_z})")
    print(l_x**2 + l_y**2 + l_z**2)

    results = {'sat_coords' : {'x_s': x_s, 'y_s': y_s, 'z_s': z_s}, 'time' : t_0, 'xi' : xi_r, 'angle' : {'l_x': l_x, 'l_y': l_y, 'l_z': l_z}}

    return t_b, results


def getAllLocationAngleCoords(t_b):
    '''
        Задача 2
        Периодичное отслеживание спутника до его выхода из зоны локации'''

    t_b = func_3(t_b)
    xi = func_1_1_reversed(t_b)
    x_s, y_s, z_s = func_1(xi, t_b)

    while func_2(x_s, y_s, z_s):
        xi_r = func_4(t_b)
        print(f"Eccentrisity anomaly =\t\t{xi_r}")
        l_x, l_y, l_z = func_5(t_b, xi_r)
        # l_x, l_y, l_z = func_5(t_b, xi)
        print(f"Station angle coords =\t\t({l_x},\t{l_y},\t{l_z})")

        t_b = func_3(t_b)
        xi = func_1_1_reversed(t_b)
        x_s, y_s, z_s = func_1(xi, t_b)


if __name__ == "__main__":

    start_time = time.time()

    t_b, _ = getFirstLocationAngleCoords()

    getAllLocationAngleCoords(t_b)

    # xi_new = func_1_1_reversed(t_b)
    # print(f"Эксцентрическая аномалия:\t{xi_new}\tот времени:\t{t_b}")

    print(f"\nProgramm works {time.time() - start_time} seconds.")
