import numpy as np
import csv
import pprint


def csv_dict_reader(file_obj: csv) -> dict:
    """ Read a CSV file using csv.DictReader and return panel_parameters"""

    reader = csv.DictReader(file_obj)

    panel_parameters = {'width_rib': list(), 'web_string': list(), 'thick_web': list(),
                        'boom_string': list(), 'thick_boom': list(),
                        'round_string': list(), 'm': list(), 'thick_cover': list(),
                        'S11': list(), 'S22': list(), 'sections': list(), 'step_string': list(),
                        'step_rib': list()}

    for line in reader:
        for k in list(panel_parameters.keys()):
            panel_parameters[k].append(line[k])

    return panel_parameters


def effective_thickness(web_string: list, thick_web: list, boom_string: list, thick_boom: list,
                        round_string: list, thick_cover: list, step_string: list) -> list:
    """" Приведеная толщина обшивки исходя из заданных параметров панели """

    EffectiveThickness = list()

    for i in range(0, len(web_string)):
        # Внешний радиус скругления стрингера
        round_external = float(round_string[i]) + (float(thick_web[i]) + float(thick_boom[i])) / 2

        # Определение площади элементов сжатого участка панели
        F_cover = (float(step_string[0]) * float(thick_cover[i]))
        F_boom = (float(boom_string[i]) - round_external) * float(thick_boom[i])
        F_web = (float(web_string[i]) - round_external) * float(thick_web[i])
        F_round = (np.pi * round_external ** 2 - np.pi * float(round_string[i]) ** 2)
        E = ((F_cover + 2 * (F_boom + F_web + F_round)) / float(step_string[0]))

        EffectiveThickness.append(E)

    return EffectiveThickness


def ratio_string(web_string: list, thick_web: list, boom_string: list, thick_boom: list,
                 round_string: list, thick_cover: list, step_string: list) -> list:
    """ Соотношение эффективной площади стрингеров к эффективной площади панели """

    ratio = list()

    for i in range(0, len(web_string)):
        # Внешний радиус скругления стрингера
        round_external = float(round_string[i]) + (float(thick_web[i]) + float(thick_boom[i])) / 2

        # Определение площади элементов сжатого участка панели
        F_cover = (float(step_string[0]) * float(thick_cover[i]))
        F_boom = (float(boom_string[i]) - round_external) * float(thick_boom[i])
        F_web = (float(web_string[i]) - round_external) * float(thick_web[i])
        F_round = (np.pi * round_external ** 2 - np.pi * float(round_string[i]) ** 2)
        r = ((2 * (F_boom + F_web + F_round)) / (F_cover + 2 * (F_boom + F_web + F_round)))

        ratio.append(r)

    return ratio


def j_panel(web_string: list, thick_web: list, boom_string: list, thick_boom: list,
            round_string: list, thick_cover: list, step_string: list) -> list:
    """ Расчет момента инерции сжатого участка панели """

    InertiaMoment = list()

    for i in range(0, len(web_string)):
        # Внешний радиус скругления стрингера
        round_external = float(round_string[i]) + (float(thick_web[i]) + float(thick_boom[i])) / 2
        # Расчет собственных момента инерции элементов сжатого участка панели
        # (пролет обшивки, полка стрингера, стенка стрингера, и т.д.)

        # Собственный момент инерции полки, кг * мм2
        Jc_boom_string = (float(boom_string[i]) - round_external) * float(thick_boom[i]) ** 3 / 12
        # Собственный момент инерции стенки, кг * мм2
        Jc_web_string = float(thick_web[i]) * (float(web_string[i]) - round_external) ** 3 / 12
        # Собственный момент инерции четверти коружности, кг * мм2
        Jc_round = np.pi * (round_external * 2) ** 4 / 256 * \
                   (1 - (float(round_string[i]) / round_external) ** 4)
        # Собственный момент инерции участка обшивки, кг * мм2
        Jc_cover = float(step_string[0]) * float(round_string[i]) ** 3 / 12

        # Определение площади элементов сжатого участка панели
        F_cover = (float(step_string[0]) * float(thick_cover[i]))
        F_boom = (float(boom_string[i]) - round_external) * float(thick_boom[i])
        F_web = (float(web_string[i]) - round_external) * float(thick_web[i])
        F_round = (np.pi * round_external ** 2 - np.pi * float(round_string[i]) ** 2)

        # Определение центра тяжести элементов сжатого участка панели
        y_cover = float(thick_cover[i]) / 2
        y_boom = (float(thick_boom[i]) / 2 + float(thick_cover[i]))
        y_web = (float(thick_cover[i]) + round_external + (float(web_string[i]) - round_external) / 2)
        y_round = (4 / 3 / np.pi * (round_external - float(round_string[i])))

        # Определение центра тяжести панели
        y = (F_cover * y_cover + (F_boom * y_boom + F_web * y_web + F_round * y_round) * 2) / \
            (F_cover + (F_boom + F_web + F_round) * 2)

        # Момент инерции сжатого участка панели, кг * мм2
        J = (Jc_cover + F_cover * (np.abs(y - y_cover)) ** 2) + \
            (Jc_boom_string + F_boom * (np.abs(y - y_boom)) ** 2) + \
            (Jc_web_string + F_web * (np.abs(y - y_web)) ** 2) + \
            (Jc_round + F_round * (np.abs(y - y_round)) ** 2)

        InertiaMoment.append(J)

    return InertiaMoment


# %% Расчет сжатой панели кессона крыла подкрепленной Т-образными стрингерами на обшую потерю устойчивости по формуле Эйлера
def sigm_eiler(InertiaMoment: list, m: list, Ex: float, step_rib: list) -> list:
    """ Критические напряжения общей потери устойчивости в секции """

    sigm = list()

    for i in range(0, len(InertiaMoment)):
        s = int(m[i]) ** 2 * np.pi ** 2 * Ex * float(InertiaMoment[i]) / float(step_rib[0]) ** 2  # кг/мм2

        sigm.append(s)

    return sigm


def safe_coeff_buckling(sigm: list, S11: list) -> list:
    SafeCoeff = list()

    for i in range(0, len(sigm)):
        SafeCoeff.append(float(sigm[i]) / float(S11[i]))

    return SafeCoeff


def bucklingCompression(C, S11: list, S22: list, step_rib: list, step_string: list,
                        boom_string: list, thick_cover: list) -> list:
    lambd_min = list()

    for k in range(0, len(thick_cover)):

        Lx = float(step_rib[0])
        # Ly = float(step_string[0]) - 2 * float(boom_string[k])
        Ly = float(thick_cover[k])

        lambd = np.zeros((5, 5))

        Nx = float(S11[k]) * float(thick_cover[k])
        Ny = float(S22[k]) * float(thick_cover[k])

        for i in range(1, 6):
            for j in range(1, 6):
                lambd[i - 1, j - 1] = np.pi ** 2 * (
                            C[0, 0] * (i / Lx) ** 4 + 2 * (C[0, 1] + 2 * C[5, 5]) * (i / Lx) ** 2 * (j / Ly) ** 2 + C[
                        1, 1] * (j / Ly) ** 4) / (Nx * (i / Lx) ** 2 + Ny * (j / Ly) ** 2)

        lambd_min.append(lambd.min())

    return lambd_min


def main_panel(panel_parameters: dict, Ex, Ey, Gxy, NUxy, C):
    EffectiveThickness = effective_thickness(web_string=panel_parameters['web_string'],
                                             thick_web=panel_parameters['thick_web'],
                                             boom_string=panel_parameters['boom_string'],
                                             thick_boom=panel_parameters['thick_boom'],
                                             round_string=panel_parameters['round_string'],
                                             thick_cover=panel_parameters['thick_cover'],
                                             step_string=panel_parameters['step_string'])

    ratio = ratio_string(web_string=panel_parameters['web_string'],
                         thick_web=panel_parameters['thick_web'],
                         boom_string=panel_parameters['boom_string'],
                         thick_boom=panel_parameters['thick_boom'],
                         round_string=panel_parameters['round_string'],
                         thick_cover=panel_parameters['thick_cover'],
                         step_string=panel_parameters['step_string'])

    InertiaMoment = j_panel(web_string=panel_parameters['web_string'],
                            thick_web=panel_parameters['thick_web'],
                            boom_string=panel_parameters['boom_string'],
                            thick_boom=panel_parameters['thick_boom'],
                            round_string=panel_parameters['round_string'],
                            thick_cover=panel_parameters['thick_cover'],
                            step_string=panel_parameters['step_string'])

    sigm = sigm_eiler(InertiaMoment=InertiaMoment,
                      m=panel_parameters['m'],
                      Ex=Ex,
                      step_rib=panel_parameters['step_rib'])

    SafeCoeff = safe_coeff_buckling(sigm=sigm, S11=panel_parameters['S11'])

    lambd_min = bucklingCompression(C=C,
                                    S11=panel_parameters['S11'],
                                    S22=panel_parameters['S22'],
                                    step_rib=panel_parameters['step_rib'],
                                    step_string=panel_parameters['step_string'],
                                    boom_string=panel_parameters['boom_string'],
                                    thick_cover=panel_parameters['thick_cover'])

    return EffectiveThickness, ratio, SafeCoeff, lambd_min

# %% Расчет сжатой панели кессона крыла подкрепленной Т-образными стрингерами на местную потерю устойчивости (расчет неподкрепленного участка обшивки)
