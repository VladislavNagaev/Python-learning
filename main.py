import numpy as np
import csv
import datetime

import stiffness_matrix as sm
import panel_analyzing as pa
import float_converter as fc

csv_path = "stiffness_matrix.csv"
with open(csv_path, "r") as f_obj:
    laminate_parameters = sm.csv_dict_reader(f_obj)
Ex, Ey, Gxy, NUxy, C = sm.effective_property_of_laminate(laminate_parameters=laminate_parameters)

csv_path2 = "panel_analyzing.csv"
with open(csv_path2, "r") as f_obj:
    panel_parameters = pa.csv_dict_reader(f_obj)
EffectiveThickness, ratio, SafeCoeff, lambd_min = pa.main_panel(panel_parameters=panel_parameters,
                                                                Ex=Ex, Ey=Ey, Gxy=Gxy,
                                                                NUxy=NUxy, C=C)

EffectiveThickness = fc.float_converter_1(convert_list=EffectiveThickness)
ratio = fc.float_converter_2(convert_list=ratio)
SafeCoeff = fc.float_converter_1(convert_list=SafeCoeff)
lambd_min = fc.float_converter_1(convert_list=lambd_min)

Ex = float("{0:.1f}".format(Ex))
Ey = float("{0:.1f}".format(Ey))
Gxy = float("{0:.1f}".format(Gxy))
NUxy = float("{0:.3f}".format(NUxy))

path = "stiffness_matrix_output.csv"
with open(path, "w", newline='') as out_file:
    writer = csv.writer(out_file)
    # writer.writerow((list(laminate_parameters.keys())))
    writer.writerow(('number, orientation, E1, E2, G, thickness, nu21, , Ex, Ey, Gxy, NUxy'.split(",")))

    n = len(laminate_parameters['number'])

    for i in range(0, n):
        if i == 0:
            writer.writerow(
                (laminate_parameters['number'][i], laminate_parameters['orientation'][i], laminate_parameters['E1'][i],
                 laminate_parameters['E2'][i], laminate_parameters['G'][i], laminate_parameters['thickness'][i],
                 laminate_parameters['nu21'][i], ' ', Ex, Ey, Gxy, NUxy))
        else:
            writer.writerow(
                (laminate_parameters['number'][i], laminate_parameters['orientation'][i], laminate_parameters['E1'][i],
                 laminate_parameters['E2'][i], laminate_parameters['G'][i], laminate_parameters['thickness'][i],
                 laminate_parameters['nu21'][i]))


path = "panel_analyzing_output.csv"
with open(path, "w", newline='') as out_file:
    writer = csv.writer(out_file)
    writer.writerow(('EffectiveThickness, ratio, SafeCoeff, lambd_min'.split(",")))

    n = len(EffectiveThickness)

    for i in range(0, n):
        writer.writerow((EffectiveThickness[i], ratio[i], SafeCoeff[i], lambd_min[i]))

# print(datetime.date.today())
# print('\nEx = {}, \nEy = {}, \nGxy = {}, \nNUxy = {}'
#       .format('%.3f' % Ex, '%.3f' % Ey, '%.3f' % Gxy, '%.3f' % NUxy))
# print('\nРеальная приведенная толщина обшивки в секции (мм) \n{}'.format(EffectiveThickness))
# print('\nСоотношение эффективной площади стрингеров к эффективной площади панели (мм) \n{}'.format(ratio))
# print('\nКоэффициент запаса общей потери устойчивости в секции \n{}'.format(SafeCoeff))
# print('\nКоэффициент запаса местной потери устойчивости в секции \n{}'.format(lambd_min))
