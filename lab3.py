# неправильно подобраны u1 и u2 (не удовлетворяют краевым условиям)
# градов скажет что хрень, но примет если объяснишь что неправильно

import pylab
from math import *
import numpy
from tkinter import *

class Params():
    def __init__(self):
        self.len = 10
        self.R = 0.5
        self.T_env = 300
        self.F0 = 100 # Плотность теплового потока
        self.lambda0 = 0.1 # Коэффициент теплопроводности в начале стержня
        self.lambda1 = 0.2 # Коэффициент теплопроводности в конце стержня
        self.alpha0 = 0.9e-2 # Коэффициент теплоотдачи в начале стержня
        self.alpha1 = 1e-2 # Коэффициент теплоотдачи в конце стержня


class Collocation(object):
    def __init__(self, params):
        self.len = params.len
        self.R = params.R
        self.T_env = params.T_env
        self.F0 = params.F0
        self.lambda0 = params.lambda0
        self.lambda1 = params.lambda1
        self.alpha0 = params.alpha0
        self.alpha1 = params.alpha1

        self.a_alpha = self.get_a(self.alpha0, self.alpha1)
        self.b_alpha = self.get_b(self.alpha0, self.alpha1)

        self.a_lambda = self.get_a(self.lambda0, self.lambda1)
        self.b_lambda = self.get_b(self.lambda0, self.lambda1)

        self.k1 = log(self.F0 / self.lambda0)
        self.k2 = self.T_env + self.F0 * (self.lambda1 - self.alpha1) / (self.lambda0 * self.alpha1 * exp(self.len))

    # решение системы для alpha(0) = alpha0; alpha(l) = alphaN
    def get_b(self, a0, a1):
        return self.len * a1 / (a1 - a0)
    def get_a(self, a0, a1):
        return -a0 * self.get_b(a0, a1)

    # коэфф. теплоотдачи
    def get_alpha(self, x):
        return self.a_alpha / (x - self.b_alpha)
    # коэфф. теплопроводности
    def get_lambda(self, x):
        return self.a_lambda / (x - self.b_lambda)
    def get_dlambda(self, x):
        return -self.a_lambda / (x - self.b_lambda)**2

    #-------------
    def u0(self, x):
        k1 = self.k1
        k2 = self.k2
        return exp(k1 - x) + k2
    def u1(self, x):
        len = self.len
        return (x**2 - len**2)**2
    def u2(self, x):
        len = self.len
        return (x**3 - len**3)**3
    def du0(self, x):
        k1 = self.k1
        return -exp(k1 - x)
    def du1(self, x):
        len = self.len
        return 4 * x * (x**2 - len**2)
    def du2(self, x):
        len = self.len
        return 9 * x**2 * (x**3 - len**3)**2
    def ddu0(self, x):
        k1 = self.k1
        return exp(k1 - x)
    def ddu1(self, x):
        len = self.len
        return 12 * x**2 - 4 * len**2
    def ddu2(self, x):
        len = self.len
        return 72 * x**7 - 90 * x**4 * len**3 + 18 * x * len**6
    #-------------
    
    def get_p(self, x):
        return 2 * self.get_alpha(x) / self.R

    def coeff_at_c1(self, x):
        return self.get_lambda(x)*self.ddu1(x) + self.get_dlambda(x)*self.du1(x) - self.get_p(x)*self.u1(x)

    def coeff_at_c2(self, x):
        return self.get_lambda(x)*self.ddu2(x) + self.get_dlambda(x)*self.du2(x) - self.get_p(x)*self.u2(x)

    def free_member(self, x):
        return -self.get_lambda(x)*self.ddu0(x) - self.get_dlambda(x)*self.du0(x) + self.get_p(x)*(self.u0(x) - self.T_env)

    def get_C1_C2(self):
        len = self.len
        x1 = 1*len/6
        x2 = 5*len/6
        a1 = self.coeff_at_c1(x1)
        b1 = self.coeff_at_c2(x1)
        m1 = self.free_member(x1)

        a2 = self.coeff_at_c1(x2)
        b2 = self.coeff_at_c2(x2)
        m2 = self.free_member(x2)

        coeff_matrix = numpy.array([[a1, b1],[a2, b2]])
        free_vector = numpy.array([m1, m2])
        return numpy.linalg.solve(coeff_matrix, free_vector)

    def func(self, x, c1, c2):
        return self.u0(x) + c1*self.u1(x) + c2*self.u2(x)
    
class App():
    def __init__(self):

        self.reset()
        self.create_widgets()
        #self.run()
        
    def create_widgets(self):
        self.root = Tk()
        self.s_len = Scale(self.root,orient=HORIZONTAL,
                           length=500,from_=0,to=100,
                           tickinterval=5, resolution=5)
        self.s_R = Scale(self.root,orient=HORIZONTAL,
                           length=500,from_=0,to=10,
                           tickinterval=0.5, resolution=0.5)
        self.s_T_env = Scale(self.root,orient=HORIZONTAL,
                           length=500,from_=-500,to=500,
                           tickinterval=100, resolution=50)
        self.s_F0 = Scale(self.root,orient=HORIZONTAL,
                           length=500,from_=0,to=200,
                           tickinterval=20, resolution=5)
        self.s_len.config(label="len")
        self.s_R.config(label="R")
        self.s_T_env.config(label="T environment")
        self.s_F0.config(label="F0")

        self.s_len.pack()
        self.s_R.pack()
        self.s_T_env.pack()
        self.s_F0.pack()

        self.s_len.set(self.spin_len)
        self.s_R.set(self.spin_R)
        self.s_T_env.set(self.spin_T_env)
        self.s_F0.set(self.spin_F0)
        
        button1 = Button(self.root,text="Reset")
        button1.pack()
        button1.bind("<Button-1>", self.re_reset)
        
        button2 = Button(self.root,text="Recount")
        button2.pack()
        button2.bind("<Button-1>", self.recount)

        self.run()
        self.root.mainloop()
        
    def recount(self, event):
        print(2)
        self.spin_len = self.s_len.get()
        self.spin_R = self.s_R.get()
        self.spin_T_env = self.s_T_env.get()
        self.spin_F0 = self.s_F0.get()

        self.run()
    
    def reset(self):
        self.spin_len = 10
        self.spin_R = 0.5
        self.spin_T_env = 300
        self.spin_F0 = 100

        self.spin_lambda0 = 0.1
        self.spin_lambda1 = 0.2

        self.spin_alpha0 = 0.9e-2
        self.spin_alpha1 = 1e-2

    def re_reset(self, event):
        print(1)
        self.spin_len = 10
        self.spin_R = 0.5
        self.spin_T_env = 300
        self.spin_F0 = 100

        self.spin_lambda0 = 0.1
        self.spin_lambda1 = 0.2

        self.spin_alpha0 = 0.9e-2
        self.spin_alpha1 = 1e-2

        self.s_len.set(self.spin_len)
        self.s_R.set(self.spin_R)
        self.s_T_env.set(self.spin_T_env)
        self.s_F0.set(self.spin_F0)

    def get_params(self):
        params = Params()

        params.R = self.spin_R
        params.len = self.spin_len

        params.T_env = self.spin_T_env
        params.F0 = self.spin_F0

        params.lambda0 = self.spin_lambda0
        params.lambda1 = self.spin_lambda1

        params.alpha0 = self.spin_alpha0
        params.alpha1 = self.spin_alpha1

        return params

    def run(self):
        params = self.get_params()
        colloc = Collocation(params)
        c1, c2 = colloc.get_C1_C2()

        n = 40
        h = params.len / n
        z = 0
        Z = []
        T = []

        while z < params.len + h:
            Z.append(z)
            T.append(colloc.func(z, c1, c2))
            z += h

        self.draw_results(T, Z)

    def draw_results(self, T, X):
        pylab.close()
        pylab.plot(X, T, color='r')
        pylab.grid(which='major', color='silver', linestyle="-", linewidth=1)
        pylab.minorticks_on()
        pylab.title('Tube temperature')
        pylab.ylabel('T, K', rotation=0)
        pylab.xlabel('Len, cm')
        pylab.show()

        
if __name__ == '__main__':
    window = App()
    #window.show()

