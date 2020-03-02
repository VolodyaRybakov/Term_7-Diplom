from math import pi, sqrt, sin, cos
import scipy.constants as constants

import get_data

INTEGRITY_CONSTANT = 0
LIGHT_SPEED = constants.c
GM = 3.98603 * (10**14)
OMEGA = 7.3 / (10**5)

data = get_data.getData('irkutsk', 'RadioAstron')

fi_0 = data['longitude']                # fi_0
fi = 0                                  # fi
r_0 = 6371302 + data['radial_distance']  # R_0
a = data['a']                           # a
e = data['eccentricity']                # e
theta_0 = data['latitude']              # theta_0
theta = data['inclination']             # theta
psi = 0                                 # psi

cos_psi = cos(psi)
sin_psi = sin(psi)
sin_theta = sin(theta)
cos_theta = cos(theta)
sin_theta_0 = sin(theta_0)
cos_theta_0 = cos(theta_0)


def __findTurningPeriod__(a: float) -> float:
    '''
        The period of rotation of an artificial earth satellite in an elliptical orbit

        parameters
        ----------
        a: float
            large semi-axis of elliptical orbit

        returns
        -------
            function from 'a': float'''

    return 2 * pi * sqrt(a**3 / GM)


def func_1(xi: float, time_xi: float) -> list([float]):
    '''
        Функция 1
        Закон движения спутника по орбите в топоцентрической системе отсчёта
        "findIndexingVectorForSatellite"'''
    def foo_x_s(xi: float, time_xi: float):
        complicated_argument = OMEGA * time_xi + fi_0 - fi
        cos_comp = cos(complicated_argument)
        sin_comp = sin(complicated_argument)
        first_braces = (cos_psi * cos_comp + sin_psi * cos_theta *
                        sin_comp) * cos_theta_0 - sin_psi * sin_theta * sin_theta_0
        second_braces = (sin_psi * cos_comp - cos_psi * cos_theta *
                         sin_comp) * cos_theta_0 + cos_psi * sin_theta * sin_theta_0
        x_s = a * first_braces * (cos(xi) - e) - a * \
            sqrt(1-(e**2)) * second_braces * sin(xi)
        return x_s

    def foo_y_s(xi: float, time_xi: float):
        complicated_argument = OMEGA * time_xi + fi_0 - fi
        cos_comp = cos(complicated_argument)
        sin_comp = sin(complicated_argument)
        first_braces = sin_psi * cos_theta * cos_comp - cos_psi * sin_comp
        second_braces = cos_psi * sin_comp + cos_psi * cos_theta * cos_comp  # sin_comp
        y_s = a * (cos(xi) - e) * first_braces + a * \
            sqrt(1 - (e**2)) * second_braces * sin(xi)
        return y_s

    def foo_z_s(xi: float, time_xi: float):
        complicated_argument = OMEGA * time_xi + fi_0 - fi
        cos_comp = cos(complicated_argument)
        sin_comp = sin(complicated_argument)
        first_braces = (sin_psi * cos_theta * sin_comp + cos_psi * cos_comp) * \
            sin_theta_0 + sin_psi * sin_theta * cos_theta_0  # sin_theta_0
        second_braces = (cos_psi * cos_theta * sin_comp - sin_psi *
                         cos_comp) * sin_theta_0 + cos_psi * sin_theta * cos_theta_0
        z_s = a * first_braces * (cos(xi) - e) + a * \
            sqrt(1 - (e**2)) * second_braces * sin(xi) - r_0
        return z_s

    x_s = foo_x_s(xi, time_xi)
    y_s = foo_y_s(xi, time_xi)
    z_s = foo_z_s(xi, time_xi)

    return (x_s, y_s, z_s)


def func_1_1(e_anomaly: float) -> float:
    '''
        Функция 1
        Нахждение времени от эксцентрической аномалии
        "timeFromEAnomaly"'''

    period = __findTurningPeriod__(a)
    return INTEGRITY_CONSTANT + (period / (2 * pi)) * (e_anomaly - sin(e_anomaly))


def func_1_1_reversed(e_anomaly, t):
    period = __findTurningPeriod__(a)
    return INTEGRITY_CONSTANT + (period / (2 * pi)) * (e_anomaly - sin(e_anomaly)) - t


def func_2(x_s: float, y_s: float, z_s: float) -> bool:
    '''
        Функция 2
        Определяет вхождение спутника в зону отслеживания
        "accuracyChecking"'''
    res = z_s - sqrt(x_s**2 + y_s**2 + z_s**2) * cos(1.22)
    return res > 0


def func_3(t_0: float):
    '''
        Функция 3
        Время посылки ладерного импульса
        "getSendingLaserImpulsTime"'''
    t_b = t_0 + 0.1 + 1
    return t_b


def func_4(e_anomaly: float, t_b: float) -> float:
    '''
        Функция 4
        Трансцендентное уравнение относительно кси_r
        "eAnomalyEquation"'''
    cos_e_anomaly = cos(e_anomaly)
    sin_e_anomaly = sin(e_anomaly)

    period = __findTurningPeriod__(a)
    # e_anomaly_b = 0

    a2R_0 = 2 * a * r_0

    complicated_argument = OMEGA * func_1_1(e_anomaly) + fi_0 - fi

    result = r_0 ** 2 - \
        -a2R_0 * \
        (cos_e_anomaly - e) * \
        (cos_theta_0 * sin_theta * sin_psi +
            sin_theta_0 *
            (cos_psi * cos(complicated_argument) +
                cos_theta * cos_psi * sin(complicated_argument)
             )
         ) + \
        -a2R_0 * sqrt(1 - (e**2)) * sin_e_anomaly * \
        (cos_theta_0 * sin_theta * cos_psi -
            sin_theta_0 *
            (sin_psi * cos(complicated_argument) -
                cos_theta * cos_psi * sin(complicated_argument)
             )
         ) + \
        a ** 2 * (1 - e * cos_e_anomaly) ** 2 - \
        -LIGHT_SPEED**2 * \
        (INTEGRITY_CONSTANT +
            ((period / (2*pi)) *
                (e_anomaly - e * sin_e_anomaly)
             ) - t_b
         )**2

    return result


def func_5(t_b: float, e_anomaly_r: float) -> list([float]):
    '''
        Функция 5
        Угловые координаты направления посылки лазерного импульса
        "findAngleCoords"'''

    delta_t = func_1_1(e_anomaly_r) - t_b
    complicated_argument = OMEGA * t_b + fi_0 - fi

    A = 1 / sqrt(
        LIGHT_SPEED**2 * delta_t**2 +
        2 * OMEGA * a * r_0 * delta_t *
        (
            (
                cos_psi * sin(complicated_argument) -
                cos_theta * sin_psi * cos(complicated_argument)
            ) * (cos(e_anomaly_r) - e) -
            sqrt(1 - e**2) *
            (
                sin_psi * sin(complicated_argument) +
                cos_theta * cos_psi * cos(complicated_argument)
            ) * sin(e_anomaly_r)
        ) * sin_theta +
        OMEGA**2 * r_0**2 * delta_t**2 * sin_theta_0**2
    )

    print(f"A = \t{A}")

    l_x = \
        A * \
        (
            a *
            (
                (
                    cos_psi * cos(complicated_argument) +
                    cos_theta * sin_psi * sin(complicated_argument)
                ) * cos_theta_0 - sin_theta * sin_theta_0 * sin_psi
            ) * (cos(e_anomaly_r) - e) -
            a * sqrt(1 - e**2) *
            (
                (
                    sin_psi * cos(complicated_argument) -
                    cos_theta * cos_psi * sin(complicated_argument)
                ) * cos_theta_0 + sin_theta * sin_theta_0 * cos_psi
            ) * sin(e_anomaly_r)
        )

    l_y = \
        (-A) * \
        (
            a *
            (
                cos_psi * sin(complicated_argument) -
                cos_theta * sin_psi * cos(complicated_argument)
            ) * (cos(e_anomaly_r) - e) -
            a * sqrt(1 - e**2) *
            (
                sin_psi * sin(complicated_argument) +
                cos_theta * cos_psi * cos(complicated_argument)
            ) * sin(e_anomaly_r) + OMEGA * r_0 * delta_t * sin_theta_0
        )

    l_z = \
        A * \
        (
            a *
            (
                (
                    cos_psi * cos(complicated_argument) +
                    cos_theta * sin_psi * sin(complicated_argument)
                ) * sin_theta_0 + sin_theta * cos_theta_0 * sin_psi
            ) * (cos(e_anomaly_r) - e) -
            a * sqrt(1 - e**2) *
            (
                (
                    sin_psi * cos(complicated_argument) -
                    cos_theta * cos_psi * sin(complicated_argument)
                ) * sin_theta_0 - sin_theta * cos_theta_0 * cos_psi
            ) * sin(e_anomaly_r) - r_0
        )

    return (l_x, l_y, l_z)