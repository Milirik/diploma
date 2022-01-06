from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


class SpaceWidget(QWidget):

	def __init__(self, parent=None):
		QWidget.__init__(self, parent)

		self.canvas = FigureCanvas(Figure())

		vertical_layout = QVBoxLayout()
		vertical_layout.addWidget(self.canvas)

		self.canvas.axes = self.canvas.figure.add_subplot(projection='3d')
		self.setLayout(vertical_layout)

		# t = np.arange(0, 10*np.pi, 0.01)
		# x, y, z = np.sin(t), np.cos(t), 0.3*t


		# Spiral=self.canvas.axes.plot(x, y, z)[0]   # просто добавляем 3 координату
		# Vertical=self.canvas.axes.plot([0, 0], [0, 0], [min(z)-0.5, max(z)+0.5],color=[0, 0, 0],linestyle='dashed')[0]

		# Point = self.canvas.axes.plot(x[0],y[0],z[0],color=[1, 0, 0],marker='o',markersize=10)[0]
		# Point2 = self.canvas.axes.plot(0,0,z[0],color=[0, 1, 0],marker='o')[0]
		# AB = self.canvas.axes.plot([0, x[0]],[0, y[0]],[z[0], z[0]],color=[0, 0, 0],linestyle='dotted')[0]

		# dt=0.01
		# def NewPoints(i):
		# 	Point.set_data_3d(x[i], y[i], z[i])   # set_data_3d. добавили 3д.
		# 	Point2.set_data_3d(0, 0, z[i])
		# 	AB.set_data_3d([0, x[i]], [0, y[i]], [z[i], z[i]])
		# 	return [Point,Point2,AB]

		# a = FuncAnimation(self.canvas.figure, NewPoints, interval=dt*1000, blit=True, frames=len(x))

		# self.canvas.draw()