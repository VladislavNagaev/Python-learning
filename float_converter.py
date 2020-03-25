def float_converter_3(convert_list: list) -> list:
    get_list = list()

    for i in range(0, len(convert_list)):
        get_list.append("{0:.3f}".format(convert_list[i]))

    return get_list


def float_converter_2(convert_list: list) -> list:
    get_list = list()

    for i in range(0, len(convert_list)):
        get_list.append("{0:.2f}".format(convert_list[i]))

    return get_list


def float_converter_1(convert_list: list) -> list:
    get_list = list()

    for i in range(0, len(convert_list)):
        get_list.append("{0:.1f}".format(convert_list[i]))

    return get_list
