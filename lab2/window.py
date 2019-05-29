from PyQt5.QtWidgets import *
from PyQt5.uic import loadUi
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

from circuit import *
from params import *


class App(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        loadUi("qt_designer.ui", self)

        self.setWindowTitle("Моделирование процессов в схеме")
        self.addToolBar(NavigationToolbar(self.MplWidget.canvas, self))

        self.btn_reset.clicked.connect(self.reset)
        self.btn_run.clicked.connect(self.run)

        self.reset()

    def reset(self):
        self.spin_R.setValue(0.35)
        self.spin_Tw.setValue(2000)

        self.spin_Ck.setValue(150)
        self.spin_Rk.setValue(0.5)
        self.spin_Lk.setValue(60)
        self.spin_Le.setValue(12)

        self.spin_Uc0.setValue(1500)
        self.spin_I0.setValue(1)

        self.spin_end.setValue(600)
        self.spin_step.setValue(1)

        self.checkBoxRp.setChecked(False)

    def get_params(self):
        """ не переводим только сантиметры """
        params = Params()

        params.R = self.spin_R.value()
        params.Tw = self.spin_Tw.value()

        params.Ck = self.spin_Ck.value() * 1e-6
        params.Rk = self.spin_Rk.value()
        params.Lk = self.spin_Lk.value() * 1e-6
        params.Le = self.spin_Le.value()

        params.Uc0 = self.spin_Uc0.value()
        params.I0 = self.spin_I0.value()

        params.end = self.spin_end.value() * 1e-6
        params.step = self.spin_step.value() * 1e-6

        params.remove_Rp = self.checkBoxRp.isChecked()

        return params

    def run(self):
        params = self.get_params()
        circuit = Circuit(params)
        time, Uc, I, Rp = circuit.solve_system()

        Up = [i*r for i, r in zip(I, Rp)]
        self.draw_results(time * 1e+6, Uc, I, Rp, Up)

    def draw_results(self, time, Uc, I, Rp, Up):
        ax1 = self.MplWidget.canvas.axes1
        ax2 = self.MplWidget.canvas.axes2
        ax3 = self.MplWidget.canvas.axes3
        ax4 = self.MplWidget.canvas.axes4

        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()

        ax1.plot(time, Uc)
        ax1.yaxis.set_label_coords(-0.055, 1.05)
        ax1.set_ylabel('Uc, В', rotation=0)
        ax1.xaxis.set_label_coords(1.05, -0.040)
        ax1.set_xlabel('t, мкс')
        ax1.grid(True)

        ax2.plot(time[1:], Rp[1:])
        ax2.yaxis.set_label_coords(-0.055, 1.05)
        ax2.set_ylabel('Rp, Ом', rotation=0)
        ax2.xaxis.set_label_coords(1.05, -0.040)
        ax2.set_xlabel('t, мкс')
        ax2.grid(True)

        ax3.plot(time, I)
        ax3.yaxis.set_label_coords(-0.055, 1.10)
        ax3.set_ylabel('I, А', rotation=0)
        ax3.xaxis.set_label_coords(1.05, -0.040)
        ax3.set_xlabel('t, мкс')
        ax3.grid(True)

        ax4.plot(time[1:], Up[1:])
        ax4.yaxis.set_label_coords(-0.055, 1.05)
        ax4.set_ylabel('Up, В', rotation=0)
        ax4.xaxis.set_label_coords(1.05, -0.040)
        ax4.set_xlabel('t, мкс')
        ax4.grid(True)

        self.MplWidget.canvas.draw()

