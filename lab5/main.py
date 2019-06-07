from data_getter import Data
from method import fixed_point_iteration
from math import fabs, exp
from numpy import arange
from copy import deepcopy
from matplotlib import pyplot as plot


def left_boundary_conditions():
    # relative to the current temperature
    X_half = Data.Xn_plus_half(Data.x0)
    p_half = Data.p(Data.x0 + Data.h / 2)
    p0 = Data.p(Data.x0)

    f0 = Data.f(Data.x0)
    f1 = Data.f(Data.x0 + Data.h)

    C0 = Data.C(Data.x0)
    C_half = Data.C(Data.x0 + Data.h / 2)

    K0 = Data.tau * (X_half / Data.h + Data.h / 8 * p_half + Data.h / 4 * p0) + \
        Data.h / 4 * C0 + Data.h / 8 * C_half

    M0 = Data.tau * (Data.h / 8 * p_half - X_half / Data.h) + Data.h / 8 * C_half

    P0 = Data.F0 * exp(-Data.gamma*t) + Data.tau * Data.h / 8 * (3 * f0 + f1) + \
        Data.h / 4 * C0 * Data.T_past[Data.x0] + \
        Data.h / 8 * C_half * (Data.T_past[Data.x0] + Data.T_past[Data.x0 + Data.h])

    return K0, M0, P0


def right_boundary_conditions():
    # relative to the current temperature
    X_half = Data.Xn_minus_half(Data.l)
    pN = Data.p(Data.l)
    p_half = Data.p(Data.l - Data.h / 2)

    fN = Data.f(Data.l)
    fN1 = Data.f(Data.l - Data.h)

    CN = Data.C(Data.l)
    C_half = Data.C(Data.l - Data.h / 2)

    KN = - Data.tau * (X_half / Data.h + Data.alphaN + Data.h / 4 * pN + Data.h / 8 * p_half) + \
        Data.h / 4 * CN + Data.h / 8 * C_half

    MN = Data.tau * (X_half / Data.h - Data.h / 8 * p_half) + Data.h / 8 * C_half

    PN = - Data.alphaN * Data.Tenv * Data.tau - Data.tau * Data.h / 8 * (3 * fN + fN1) + \
        Data.h / 4 * CN * Data.T_past[Data.l] + \
        Data.h / 8 * C_half * (Data.T_past[Data.l] + Data.T_past[Data.l - Data.h])

    return KN, MN, PN


def calc_coefficients():
    # relative to the current temperature
    A = []
    B = []
    D = []
    F = []

    for i in arange(Data.x0, Data.l, Data.h):
        i = round(i, Data.rounding)
        An = Data.tau * Data.Xn_minus_half(i) / Data.h
        Dn = Data.tau * Data.Xn_plus_half(i) / Data.h
        Bn = An + Dn + Data.p(i) * Data.h * Data.tau + Data.C(i) * Data.h
        Fn = Data.C(i) * Data.h * Data.T_past[i] + Data.f(i) * Data.h * Data.tau

        A.append(An)
        B.append(Bn)
        D.append(Dn)
        F.append(Fn)

    return A, B, D, F


def calc_changes(a, b):
    lb = len(b)
    la = len(a)

    if lb > la:
        a, b = b, a

    diff = []
    for i in range(la):
        if lb > i:
            diff.append(fabs(b[i] - a[i]))
        else:
            diff.append(a[i])
    return diff


def get_values_from_dict(x: dict):
    return [x[item] for item in x.keys()]


def distribute_temperature(T):
    res = {}
    x = Data.x0
    for t in T:
        res.update({round(x, Data.rounding): t})
        x += Data.h

    return res


if __name__ == "__main__":
    max_Temp = 0
    T_old = get_values_from_dict(Data.T_curr)
    xs = list(Data.T_curr.keys())
    plot.plot(xs, T_old)

    #print(T_old)
    t = 0
    step = 0
    max_Tmp = None
    max_tmp = None
    time = []
    time.append(0)
    TM = [max(T_old)]
    claa = []
    F0 = [0]
    while True:
        T_new = deepcopy(T_old)
        Data.T_past = distribute_temperature(T_old)
        while True:
            T_old = deepcopy(T_new)
            time.append(time[-1] + Data.tau)
            #print(T_old)
            a, b, c, d = calc_coefficients()
            k0, m0, p0 = left_boundary_conditions()
            kN, mN, pN = right_boundary_conditions()
            T_new = fixed_point_iteration(a, b, c, d, k0, m0, p0, kN, mN, pN)
            TM.append(max(T_new))
            F0.append(Data.F0 * time[-1] * exp(-Data.gamma * time[-1]))

            diff = calc_changes(T_old, T_new)
            m = max(diff)
            if m / T_new[diff.index(m)] < Data.eps:
                Data.T_curr = distribute_temperature(T_new)
                step += 1
                # print(T_new)
                if step // 25 == 0:  # 42
                    #plot.subplot(311)
                    plot.plot(xs, T_new, "skyblue")
                    flag_max = True
                    a_index = 0

                    if max(T_new) > max_Temp:
                        max_Temp = max(T_new)
                        max_Tmp = T_new
                    max_tmp = T_new
                    #if First_red is not None:
                    #    plot.plot(xs[:len(First_red)], First_red, "g")
                    #plot.plot(xs, Second_green, "g")
                break

        diff = calc_changes(get_values_from_dict(Data.T_past), T_new)
        t += 1
        if t == 23:
            break
        m = max(diff)
        #print(m/ T_new[diff.index(m)], m/ T_new[diff.index(m)] < Data.eps)
        if m / T_new[diff.index(m)] < Data.eps:
            Data.T_past = distribute_temperature(T_old)
            break
    #print(step)
    #plot.subplot(311)

    plot.plot(xs, max_Tmp, "r")
    plot.plot(xs, max_tmp, "b")
    plot.grid()


    plot.show()

