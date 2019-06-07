from numpy import arange
from method import two_d_interpolation
import ast
import configparser

class Data:
    gamma = 0.4
    theta = 300
    #coeff = -0.2399 #-0.3964
    c0 = 0.163 #0.8
    x0 = 0
    l = 10           # Length of the rod (cm)
    R = 0.5          # Radius of the rod (cm)
    Tenv = 300       # Ambient temperature (K)
    F0 = 100         # Heat flux density (W / (cm^2 * K))
    alpha0 = 1e-2    # Heat transfer coefficient at the beginning of the rod (W / (cm^2 * K))
    alphaN = 9e-3  # Heat transfer coefficient at the end of the rod (W / (cm^2 * K))
    h = 1e-2
    tau = 0.25
    eps = 1e-2
    rounding = 4
    b_alpha = (alphaN * l) / (alphaN - alpha0)
    a_alpha = - alpha0 * b_alpha

    T_curr = {round(x, 4): 300 for x in arange(x0, l+h, h)}
    t_curr = 0

    T_past = None

    @staticmethod
    def get_table_of_heat_cap():
        config = configparser.ConfigParser()
        config.read("./config.ini")

        # Читаем некоторые значения из конфиг. файла.
        font = config.get("table_of_heat_cap", "value")
        return ast.literal_eval(font)

    @staticmethod
    def get_table_of_thermal_cond():
        config = configparser.ConfigParser()
        config.read("./config.ini")

        # Читаем некоторые значения из конфиг. файла.
        font = config.get("table_of_thermal_cond", "value")
        return ast.literal_eval(font)

    @staticmethod
    def alpha(x):
        return Data.a_alpha / (x - Data.b_alpha)

    @staticmethod
    def k(x):
        if x in Data.T_curr.keys():
            T = Data.T_curr[x]
        else:
            T = two_d_interpolation(Data.T_curr, x)
        table = Data.get_table_of_thermal_cond()
        if T in table.keys():
            return table[T]
        else:
            # interpolation
            return two_d_interpolation(table, T)

    @staticmethod
    def C(x):
        '''
        print(Data.c0 * pow(x, Data.coeff), pow(Data.theta, Data.coeff))
        return Data.c0 * pow(x, Data.coeff)/pow(Data.theta, Data.coeff)
        '''
        if x in Data.T_curr.keys():
            T = Data.T_curr[x]
        else:
            T = two_d_interpolation(Data.T_curr, x)
        table = Data.get_table_of_heat_cap()
        if T in table.keys():
            return table[T]
        else:
            # interpolation
            return two_d_interpolation(table, T)

    @staticmethod
    def Xn_minus_half(x):
        return (Data.k(x - Data.h) + Data.k(x)) / 2

    @staticmethod
    def Xn_plus_half(x):
        return (Data.k(x) + Data.k(x + Data.h)) / 2

    @staticmethod
    def p(x):
        return 2 * Data.alpha(x) / Data.R

    @staticmethod
    def f(x):
        return 2 * Data.alpha(x) / Data.R * Data.Tenv


if __name__ == "__main__":
    print(Data.C(0.49585855))


