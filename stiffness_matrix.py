import csv
import numpy as np


def csv_dict_reader(file_obj: csv) -> dict:
    """ Read a CSV file using csv.DictReader and return laminate_parameters"""

    reader = csv.DictReader(file_obj)

    laminate_parameters = {'number': list(), 'orientation': list(), 'E1': list(), 'E2': list(),
                           'G': list(), 'thickness': list(), 'nu21': list()}

    for line in reader:
        for k in list(laminate_parameters.keys()):
            laminate_parameters[k].append(line[k])

    return laminate_parameters


def nu_12(nu21: list, E1: list, E2: list) -> list:
    """Вычисление nu21"""

    n = len(nu21)
    nu12 = list()

    for i in range(0, n):
        nu12.append(float(nu21[i]) * float(E2[i]) / float(E1[i]))

    return nu12


def stiffness_matrix_of_lamina(E1: list, E2: list, G: list, nu12: list, nu21: list, orientation: list):
    """Вычисление матрицы жесткости монослоя в системе координат пакета"""

    n = len(nu21)
    Q = np.zeros((n, 3, 3))
    Qbar = np.zeros((n, 3, 3))

    for i in range(0, n):

        Qt = np.array([[float(E1[i]) / (1 - float(nu12[i]) * float(nu21[i])),
                        float(nu12[i]) * float(E1[i]) / (1 - float(nu12[i]) * float(nu21[i])), 0],
                       [float(nu12[i]) * float(E1[i]) / (1 - float(nu12[i]) * float(nu21[i])),
                        float(E2[i]) / (1 - float(nu12[i]) * float(nu21[i])), 0],
                       [0, 0, float(G[i])]])

        Q[i, :, :] = Qt

        c = np.cos(float(orientation[i]) * np.pi / 180)
        s = np.sin(float(orientation[i]) * np.pi / 180)

        T = np.array([[c ** 2, s ** 2, 2 * c * s],
                      [s ** 2, c ** 2, -2 * c * s],
                      [-c * s, c * s, c ** 2 - s ** 2]])

        Tinv = np.array([[c ** 2, s ** 2, -1 * c * s],
                         [s ** 2, c ** 2, 1 * c * s],
                         [2 * c * s, -2 * c * s, c ** 2 - s ** 2]])

        Qtt_temp = np.matmul(T, Qt)
        Qtt = np.matmul(Qtt_temp, Tinv)
        Qbar[i, :, :] = Qtt

    return Qbar


def coordinate_z(thickness: list):
    """Координаты слоев"""

    n = len(thickness)
    z = np.zeros((n + 1))

    for i in range(0, n):
        z[i + 1] = z[i] + float(thickness[i])

    z = z - np.max(z) / 2

    return z


def abd_matrix(Qbar, z):
    """Матрица жесткости"""

    n = len(z) - 1
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))

    for i in range(0, n):
        A = A + Qbar[i, :, :] * (z[i + 1] - z[i])
        B = B + 0.5 * Qbar[i, :, :] * (z[i + 1] ** 2 - z[i] ** 2)
        D = D + 1 / 3 * Qbar[i, :, :] * (z[i + 1] ** 3 - z[i] ** 3)

    C = np.zeros((6, 6))
    C = np.array((np.r_[np.c_[A, B], np.c_[B, D]]))

    for i in range(0, 6):
        for j in range(0, 6):
            if abs(C[i, j]) <= 10 ** (-10):
                C[i, j] = 0

    return C


def effective_modules(C, z):
    """Вычисление эффективных модулей"""

    S = np.linalg.inv(C)

    for i in range(0, 6):
        for j in range(0, 6):
            if abs(S[i, j]) <= 10 ** (-10):
                S[i, j] = 0

    Ex = 1 / S[0, 0] / (np.max(z) - np.min(z))
    Ey = 1 / S[1, 1] / (np.max(z) - np.min(z))
    Gxy = 1 / S[2, 2] / (np.max(z) - np.min(z))
    NUxy = np.abs(S[0, 1] / S[0, 0])

    return Ex, Ey, Gxy, NUxy, S


def effective_property_of_laminate(laminate_parameters: dict):

    nu12 = nu_12(nu21=laminate_parameters['nu21'], E1=laminate_parameters['E1'], E2=laminate_parameters['E2'])
    Qbar = stiffness_matrix_of_lamina(E1=laminate_parameters['E1'], E2=laminate_parameters['E2'],
                                      G=laminate_parameters['G'], nu12=nu12, nu21=laminate_parameters['nu21'],
                                      orientation=laminate_parameters['orientation'])
    z = coordinate_z(thickness=laminate_parameters['thickness'])
    C = abd_matrix(Qbar=Qbar, z=z)
    Ex, Ey, Gxy, NUxy, S = effective_modules(C=C, z=z)

    return Ex, Ey, Gxy, NUxy, C


if __name__ == "__main__":

    csv_path = "stiffness_matrix.csv"

    with open(csv_path, "r") as f_obj:
        laminate_parameters = csv_dict_reader(f_obj)

    Ex, Ey, Gxy, NUxy, C = effective_property_of_laminate(laminate_parameters=laminate_parameters)

    print('Ex = {}, \nEy = {}, \nGxy = {}, \nNUxy = {}'.format('%.3f' % Ex, '%.3f' % Ey, '%.3f' % Gxy, '%.3f' % NUxy))
