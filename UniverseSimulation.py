from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import FormOfUniverseSimulation
import spacewidget
import math
import numpy as np
import random as rd
from matplotlib.animation import FuncAnimation
import sympy as sp
import pprint
import time
import scipy.io as io
import pickle
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure


class PlanetSystem():
    def __init__(self, planets, spaceShip='NotEnoughGold'):
        self.planets = planets
        if spaceShip != 'NotEnoughGold':
            self.spaceShip = spaceShip

    def AddNewPlanet(self, planet):
        self.planets.append(planet)

    def AddSpaceShip(self, spaceShip):
        self.spaceShip = spaceShip

    def ReplaceSystem(self, X, Y, Z, VX, VY, VZ, X_Sh=0, Y_Sh=0, Z_Sh=0, VX_Sh=0, VY_Sh=0, VZ_Sh=0, Phi_Sh=0, F_curr=0):
        for planet, x, y, z, vx, vy, vz in zip(self.planets, X, Y, Z, VX, VY, VZ):
            planet.Replace(x, y, z, vx, vy, vz)
            planet.ReDraw()
        if (self.spaceShip):
            self.spaceShip.Replace(X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, Phi_Sh, F_curr)
            self.spaceShip.ReDraw()

    def Draw(self, axes):
        for planet in self.planets:
            planet.Draw(axes)
        if (self.spaceShip):
            self.spaceShip.Draw(axes)

    def GetMoveEquations(self):
        n = len(self.planets)
        _strX = ''
        _strY = ''
        _strZ = ''
        _strVx = ''
        _strVy = ''
        _strVz = ''
        for i in range(n):
            _strX += f'x{i}, '
            _strY += f'y{i}, '
            _strZ += f'z{i}, '
            _strVx += f'Vx{i}, '
            _strVy += f'Vy{i}, '
            _strVz += f'Vz{i}, '

        X = sp.symbols(_strX)
        Y = sp.symbols(_strY)
        Z = sp.symbols(_strZ)

        VX = sp.symbols(_strVx)
        VY = sp.symbols(_strVy)
        VZ = sp.symbols(_strVz)


        DX = [Vx for Vx in VX]
        DY = [Vy for Vy in VY]
        DZ = [Vz for Vz in VZ]

        DVX = [
            sum([
                (planet.m * (x - cur_x)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2 + (z - cur_z)**2) ** 3)
                for x, y, z, planet in zip(X, Y, Z, self.planets)
                if (x != cur_x)
            ])
            for cur_x, cur_y, cur_z, current_planet in zip(X, Y, Z, self.planets)
        ]

        DVY = [
            sum([
                (planet.m * (y - cur_y)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2 + (z - cur_z)**2) ** 3)
                for x, y, z, planet in zip(X, Y, Z, self.planets)
                if (x != cur_x)
            ])
            for cur_x, cur_y, cur_z, current_planet in zip(X, Y, Z, self.planets)
        ]

        DVZ = [
            sum([
                (planet.m * (z - cur_z)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2 + (z - cur_z)**2) ** 3)
                for x, y, z, planet in zip(X, Y, Z, self.planets)
                if (x != cur_x)
            ])
            for cur_x, cur_y, cur_z, current_planet in zip(X, Y, Z, self.planets)
        ]

        self.SpaceBodyMoveEquations = sp.lambdify([X, Y, Z, VX, VY, VZ], [DX, DY, DZ, DVX, DVY, DVZ])

        if (self.spaceShip):
            X_Sh = sp.symbols('x_Sh')
            Y_Sh = sp.symbols('y_Sh')
            Z_Sh = sp.symbols('z_Sh')

            VX_Sh = sp.symbols('Vx_Sh')
            VY_Sh = sp.symbols('Vy_Sh')
            VZ_Sh = sp.symbols('Vz_Sh')

            F_dv = sp.symbols('f_dv')
            Alpha = sp.symbols('alpha')
            Z_boost = sp.symbols('z_boost')


            DX_Sh = VX_Sh
            DY_Sh = VY_Sh
            DZ_Sh = VZ_Sh

            DVX_Sh = sum([
                (planet.m * (x - X_Sh)) / (sp.sqrt((x - X_Sh) ** 2 + (y - Y_Sh) ** 2 + (z - Z_Sh) ** 2) ** 3)
                for x, y, z, planet in zip(X, Y, Z, self.planets)
            ]) + F_dv / self.spaceShip.m * sp.cos(Z_boost)* sp.cos(Alpha)

            DVY_Sh = sum([
                (planet.m * (y - Y_Sh)) / (sp.sqrt((x - X_Sh) ** 2 + (y - Y_Sh) ** 2 + (z - Z_Sh) ** 2) ** 3)
                for x, y, z, planet in zip(X, Y, Z, self.planets)
            ]) + F_dv / self.spaceShip.m * sp.cos(Z_boost)* sp.sin(Alpha)


            DVZ_Sh = sum([
                (planet.m * (z - Z_Sh)) / (sp.sqrt((x - X_Sh) ** 2 + (y - Y_Sh) ** 2 + (z - Z_Sh) ** 2) ** 3)
                for x, y, z, planet in zip(X, Y, Z, self.planets)
            ]) + F_dv / self.spaceShip.m * sp.sin(Z_boost)


        self.SpaceShipMoveEquations = sp.lambdify([X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, X, Y, Z, VX, VY, VZ, F_dv, Alpha, Z_boost],
                                                  [DX_Sh, DY_Sh, DZ_Sh, DVX_Sh, DVY_Sh, DVZ_Sh])

    def GetStateVectors(self):
        X = np.zeros(len(self.planets))
        Y = np.zeros(len(self.planets))
        Z = np.zeros(len(self.planets))
        VX = np.zeros(len(self.planets))
        VY = np.zeros(len(self.planets))
        VZ = np.zeros(len(self.planets))

        for i in range(len(self.planets)):
            print('[0]', self.planets[i].x, self.planets[i].y, self.planets[i].z, self.planets[i].Vx, self.planets[i].Vy, self.planets[i].Vz, type(self.planets[i].x))
            X[i] = self.planets[i].x
            Y[i] = self.planets[i].y
            Z[i] = self.planets[i].z

            VX[i] = self.planets[i].Vx
            VY[i] = self.planets[i].Vy
            VZ[i] = self.planets[i].Vz


        return X, Y, Z, VX, VY, VZ

class Planet():
    def __init__(self, x0, y0, z0, Vx0, Vy0, Vz0, m, R, color):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.Vx0 = Vx0
        self.Vy0 =  Vy0
        self.Vz0 = Vz0
        self.m = m
        self.R = R
        self.color = color

        self.x = x0
        self.y = y0
        self.z = z0
        self.Vx = Vx0
        self.Vy = Vy0
        self.Vz = Vz0

        phi = np.linspace(0, 6.28, 20)
        # self.PlanetX = self.R * np.sin(phi)
        # self.PlanetY = self.R * np.cos(phi)
        # self.PlanetZ = self.R * np.

        self.TraceX = np.array([self.x])
        self.TraceY = np.array([self.y])
        self.TraceZ = np.array([self.z])


    def Replace(self, x, y, z, vx, vy, vz):
        self.x = x
        self.y = y
        self.z = z
        self.Vx = vx
        self.Vy = vy
        self.Vz = vz

        self.TraceX = np.append(self.TraceX, x)
        self.TraceY = np.append(self.TraceY, y)
        self.TraceZ = np.append(self.TraceZ, z)

    def Draw(self, axes):
        # print('x y z', self.x, self.y, self.z)
        # print('Vx Vy Vz', self.Vx, self.Vy, self.Vz)
        # print('planet R', self.R)
        self.DrawedPlanet = axes.plot(self.x, self.y, self.z, color='red', marker='o', markersize=self.R)[0]
        self.DrawedTrace = axes.plot(self.TraceX, self.TraceY, self.TraceZ, ':')[0]

    def ReDraw(self):
        self.DrawedPlanet.set_data_3d(self.x, self.y, self.z)
        self.DrawedTrace.set_data_3d(self.TraceX, self.TraceY, self.TraceZ)

class SpaceShip():
    def __init__(self, x0, y0, z0, Vx0, Vy0, Vz0, m, R, color, FlameColor, Phi, F_max):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.Vx0 = Vx0
        self.Vy0 = Vy0
        self.Vz0 = Vz0
        self.m = m
        self.R = R
        self.color = color
        self.flame_color = FlameColor
        self.Phi = Phi

        self.x = x0
        self.y = y0
        self.z = z0
        self.Vx = Vx0
        self.Vy = Vy0
        self.Vz = Vz0

        self.phi = 0
        self.F_dv = F_max
        self.F_curr = 0

        self.SpaceShipX = self.x 
        self.SpaceShipY = self.y
        self.SpaceShipZ = self.z

        self.SpaceShipFlameX = self.R * np.array([0, -3.5,-3,-4,-3,-3.5, 0])
        self.SpaceShipFlameY = self.R * np.array([0.2, 0.15, 0.1, 0, -0.1, -0.15, -0.2])

        self.TraceX = np.array([self.x])
        self.TraceY = np.array([self.y])
        self.TraceZ = np.array([self.z])

    def Replace(self, x, y, z, vx, vy, vz, phi, F_curr):
        self.x = x
        self.y = y
        self.z = z
        self.Vx = vx
        self.Vy = vy
        self.Vz = vz
        self.phi = phi
        self.F_curr = F_curr

        self.TraceX = np.append(self.TraceX, x)
        self.TraceY = np.append(self.TraceY, y)
        self.TraceZ = np.append(self.TraceZ, z)

    def Draw(self, axes):
        print('spaceShip R', self.R)

        self.DrawedSpaceShip = axes.plot(self.x, self.y, self.z, color='black', marker='o', markersize=self.R)[0]
        # self.DrawedSpaceShipFlame = axes.plot(self.x + self.SpaceShipFlameX, self.y + self.SpaceShipFlameY, self.z, color=self.flame_color)[0]
        self.DrawedTrace = axes.plot(self.TraceX, self.TraceY, self.TraceZ, ':')[0]

    def ReDraw(self):
        # RotSpaceShipX, RotSpaceShipY, RotSpaceShipZ = Rot3D(self.x, self.y, self.z, self.phi)
        self.DrawedSpaceShip.set_data_3d(self.x, self.y, self.z)

        # PFlameX = self.F_curr/self.F_dv*self.SpaceShipFlameX + self.SpaceShipX[int(len(self.SpaceShipX)/2)]
        # RotSpaceShipFlameX, RotSpaceShipFlameY, RotSpaceShipFlameZ= Rot3D(PFlameX, self.SpaceShipFlameY, self.phi)
        # self.DrawedSpaceShipFlame.set_data_3d(self.x + RotSpaceShipFlameX, self.y + RotSpaceShipFlameY, self.z + RotSpaceShipFlameZ)
        self.DrawedTrace.set_data_3d(self.TraceX, self.TraceY, self.TraceZ)

def Rot3D(X, Y, Z,phi):
    rotateX = np.array([[1, 0, 0, 0], [0, np.cos(phi), np.sin(phi), 0], [0, -np.sin(phi), np.cos(phi), 0], [0, 0, 0, 1]])
    rotateY = np.array([[np.cos(phi), 0, np.sin(phi), 0], [0, 1, 0, 0], [-np.sin(phi), 0, np.cos(phi), 0], [0, 0, 0, 1]])
    rotateZ = np.array([[np.cos(phi), -np.sin(phi), 0, 0], [np.sin(phi), np.cos(phi), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    rotate = rotateX * rotateY * rotateZ

    tmpArray = np.array([X, Y, Z, 1])
    tmp2 = rotate * tmpArray
    RotX = tmp2[0][0]
    RotY = tmp2[1][0]
    RotZ = tmp2[2][0]

    return RotX, RotY, RotZ

def DrawTheSpace(axes):
    axes.fill([-100, 100, 100, - 100], [-100, - 100, 100, 100], 'black')
    nstars = 5000
    xstars = 200 * np.random.random(nstars) - 100
    ystars = 200 * np.random.random(nstars) - 100
    brstars = 0.5 + np.random.random(nstars) / 2
    sizestars = np.random.random(nstars)
    for i in range(nstars):
        axes.plot(xstars[i], ystars[i], xstars[i], marker='o', markersize=sizestars[i], color=[brstars[i], brstars[i], brstars[i]])

global t

class SpaceWidget(QMainWindow, FormOfUniverseSimulation.Ui_MainWindow):
    def __init__(self):
        QMainWindow.__init__(self)

        self.setupUi(self)

        self.setWindowTitle("В одной далёкой галактике")

        self.StartButton.clicked.connect(self.HereAreWeGo)
        self.addToolBar(NavigationToolbar(self.SpWidget.canvas, self))



    def HereAreWeGo(self):
        global t, dt, plSystem, X, Y, Z, VX, VY, VZ, Dx, Dy, Dz, DVx, DVy, DVz, X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, Dx_Sh, Dy_Sh, Dz_Sh, DVx_Sh, DVy_Sh, DVZ_Sh, F_max, F_dv, Alpha
        
        #     Параметры массы
        FileName = self.FileName_field.text() + '.universe'
        dt = float(self.TStep_field.text())
        K = float(self.K_field.text())

        with open(FileName, 'rb') as config_dictionary_file:
            # Step 3
            config_dictionary = pickle.load(config_dictionary_file)

        plSystem = config_dictionary['PS']
        plSystem.GetMoveEquations()


        X, Y, Z, VX, VY, VZ = plSystem.GetStateVectors()

        if (plSystem.spaceShip):
            X_Sh = plSystem.spaceShip.x
            Y_Sh = plSystem.spaceShip.y
            Z_Sh = plSystem.spaceShip.z

            VX_Sh = plSystem.spaceShip.Vx
            VY_Sh = plSystem.spaceShip.Vy
            VZ_Sh = plSystem.spaceShip.Vz

            F_max=plSystem.spaceShip.F_dv
            F_dv = self.F_Bar.value()*F_max/100
            Alpha = self.Angle_Bar.value()/360*6.28+1.57
        else:
            X_Sh = 0
            Y_Sh = 0
            Z_Sh = 0
            VX_Sh = 0
            VY_Sh = 0
            VZ_Sh = 0

            F_max = 0
            F_dv = 0
            Alpha = 0


        self.SpWidget.canvas.axes.clear()
        #self.SpWidget.canvas.axes.grid(True)
        #self.SpWidget.canvas.axes.axis('scaled')
        Side = K*20

        self.SpWidget.canvas.axes.set(xlim=[-2*Side, 2*Side], ylim=[-Side, Side], zlim=[-Side, Side])
        self.SpWidget.canvas.axes.set_title('Это космос')
        self.SpWidget.canvas.axes.set_xlabel('X')
        self.SpWidget.canvas.axes.set_ylabel('Y')
        self.SpWidget.canvas.axes.set_zlabel('Z')



        t = 0.0
        plSystem.Draw(self.SpWidget.canvas.axes)
        self.SpWidget.canvas.show()

        def NewPoints(i):
            global t, dt, plSystem, X, Y, Z, VX, VY, VZ, Dx, Dy, Dz, DVx, DVy, DVz, X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, Dx_Sh, Dy_Sh, Dz_Sh, DVx_Sh, DVy_Sh, DVz_Sh, F_max, F_dv, Alpha
            t += dt
            F_dv = self.F_Bar.value()*F_max/100
            Alpha = -self.Angle_Bar.value()/360*6.28+1.57
            Z_boost = self.Z_turbo.value()/360*6.28

            print(Z_boost)
            # Методом Рунге - Кутты
            Dx1, Dy1, Dz1, DVx1, DVy1, DVz1 = plSystem.SpaceBodyMoveEquations(X, Y, Z, VX, VY, VZ)
            Dx1_Sh, Dy1_Sh, Dz1_Sh, DVx1_Sh, DVy1_Sh, DVz1_Sh = plSystem.SpaceShipMoveEquations(X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, X, Y, Z, VX, VY, VZ,
                                                                               F_dv,
                                                                               Alpha,
                                                                               Z_boost
                                                                               )
            Dx1 = np.array(Dx1)
            Dy1 = np.array(Dy1)
            Dz1 = np.array(Dz1)
            DVx1 = np.array(DVx1)
            DVy1 = np.array(DVy1)
            DVz1 = np.array(DVz1)

            Dx1_Sh = np.array(Dx1_Sh)
            Dy1_Sh = np.array(Dy1_Sh)
            Dz1_Sh = np.array(Dz1_Sh)

            DVx1_Sh = np.array(DVx1_Sh)
            DVy1_Sh = np.array(DVy1_Sh)
            DVz1_Sh = np.array(DVz1_Sh)


            Dx2, Dy2, Dz2, DVx2, DVy2, DVz2 = plSystem.SpaceBodyMoveEquations(X + Dx1 / 2 * dt, Y + Dy1 / 2 * dt, Z + Dz1 / 2 * dt,
                                                                   VX + DVx1 / 2 * dt,
                                                                   VY + DVy1 / 2 * dt,
                                                                   VZ + DVz1 / 2 * dt,
                                                                   )

            Dx2_Sh, Dy2_Sh, Dz2_Sh, DVx2_Sh, DVy2_Sh, DVz2_Sh = plSystem.SpaceShipMoveEquations(
                X_Sh + Dx1_Sh / 2 * dt, Y_Sh + Dy1_Sh / 2 * dt, Z_Sh + Dz1_Sh / 2 * dt, VX_Sh + DVx1_Sh / 2 * dt, VY_Sh + DVy1_Sh / 2 * dt, VZ_Sh + DVz1_Sh / 2 * dt,
                X + Dx1 / 2 * dt, Y + Dy1 / 2 * dt, Z + Dz1 / 2 * dt, VX + DVx1 / 2 * dt, VY + DVy1 / 2 * dt, VZ + DVz1 / 2 * dt, F_dv, Alpha, Z_boost)
            

            Dx2 = np.array(Dx2)
            Dy2 = np.array(Dy2)
            Dz2 = np.array(Dz2)
            DVx2 = np.array(DVx2)
            DVy2 = np.array(DVy2)
            DVz2 = np.array(DVz2)
            Dx2_Sh = np.array(Dx2_Sh)
            Dy2_Sh = np.array(Dy2_Sh)
            Dz2_Sh = np.array(Dz2_Sh)
            DVx2_Sh = np.array(DVx2_Sh)
            DVy2_Sh = np.array(DVy2_Sh)
            DVz2_Sh = np.array(DVz2_Sh)

            Dx3, Dy3, Dz3, DVx3, DVy3, DVz3 = plSystem.SpaceBodyMoveEquations(X + Dx2 / 2 * dt, Y + Dy2 / 2 * dt, Z + Dz2 / 2 * dt,
                                                                   VX + DVx2 / 2 * dt,
                                                                   VY + DVy2 / 2 * dt,
                                                                   VZ + DVz2 / 2 * dt
                                                                   )

            Dx3_Sh, Dy3_Sh, Dz3_Sh, DVx3_Sh, DVy3_Sh, DVz3_Sh = plSystem.SpaceShipMoveEquations(
                X_Sh + Dx2_Sh / 2 * dt, Y_Sh + Dy2_Sh / 2 * dt, Z_Sh + Dz2_Sh / 2 * dt, VX_Sh + DVx2_Sh / 2 * dt, VY_Sh + DVy2_Sh / 2 * dt, VZ_Sh + DVz2_Sh / 2 * dt,
                X + Dx2 / 2 * dt, Y + Dy2 / 2 * dt, Z + Dz2 / 2 * dt, VX + DVx2 / 2 * dt, VY + DVy2 / 2 * dt, VZ + DVz2 / 2 * dt, F_dv, Alpha, Z_boost)

            Dx3 = np.array(Dx3)
            Dy3 = np.array(Dy3)
            Dz3 = np.array(Dz3)
            DVx3 = np.array(DVx3)
            DVy3 = np.array(DVy3)
            DVz3 = np.array(DVz3)
            Dx3_Sh = np.array(Dx3_Sh)
            Dy3_Sh = np.array(Dy3_Sh)
            Dz3_Sh = np.array(Dz3_Sh)
            DVx3_Sh = np.array(DVx3_Sh)
            DVy3_Sh = np.array(DVy3_Sh)
            DVz3_Sh = np.array(DVz3_Sh)

            Dx4, Dy4, Dz4, DVx4, DVy4, DVz4 = plSystem.SpaceBodyMoveEquations(X + Dx3 * dt, Y + Dy3 * dt, Z + Dz3 * dt, VX + DVx3 * dt,
                                                                   VY + DVy3 * dt, VZ + DVz3 * dt)

            Dx4_Sh, Dy4_Sh, Dz4_Sh, DVx4_Sh, DVy4_Sh, DVz4_Sh = plSystem.SpaceShipMoveEquations(
                X_Sh + Dx3_Sh * dt, Y_Sh + Dy3_Sh * dt, Z_Sh + Dz3_Sh * dt, VX_Sh + DVx3_Sh * dt, VY_Sh + DVy3_Sh * dt, VZ_Sh + DVz3_Sh * dt,
                X + Dx3 * dt, Y + Dy3 * dt, Z + Dz3 * dt, VX + DVx3 * dt, VY + DVy3 * dt, VZ + DVz3 * dt, F_dv, Alpha, Z_boost)

            Dx4 = np.array(Dx4)
            Dy4 = np.array(Dy4)
            Dz4 = np.array(Dz4)
            DVx4 = np.array(DVx4)
            DVy4 = np.array(DVy4)
            DVz4 = np.array(DVz4)
            Dx4_Sh = np.array(Dx4_Sh)
            Dy4_Sh = np.array(Dy4_Sh)
            Dz4_Sh = np.array(Dz4_Sh)
            DVx4_Sh = np.array(DVx4_Sh)
            DVy4_Sh = np.array(DVy4_Sh)
            DVz4_Sh = np.array(DVz4_Sh)

            X = X + dt / 6 * (Dx1 + 2 * Dx2 + 2 * Dx3 + Dx4)
            Y = Y + dt / 6 * (Dy1 + 2 * Dy2 + 2 * Dy3 + Dy4)
            Z = Z + dt / 6 * (Dz1 + 2 * Dz2 + 2 * Dz3 + Dz4)

            VX = VX + dt / 6 * (DVx1 + 2 * DVx2 + 2 * DVx3 + DVx4)
            VY = VY + dt / 6 * (DVy1 + 2 * DVy2 + 2 * DVy3 + DVy4)
            VZ = VZ + dt / 6 * (DVz1 + 2 * DVz2 + 2 * DVz3 + DVz4)

            X_Sh = X_Sh + dt / 6 * (Dx1_Sh + 2 * Dx2_Sh + 2 * Dx3_Sh + Dx4_Sh)
            Y_Sh = Y_Sh + dt / 6 * (Dy1_Sh + 2 * Dy2_Sh + 2 * Dy3_Sh + Dy4_Sh)
            Z_Sh = Z_Sh + dt / 6 * (Dz1_Sh + 2 * Dz2_Sh + 2 * Dz3_Sh + Dz4_Sh)

            VX_Sh = VX_Sh + dt / 6 * (DVx1_Sh + 2 * DVx2_Sh + 2 * DVx3_Sh + DVx4_Sh)
            VY_Sh = VY_Sh + dt / 6 * (DVy1_Sh + 2 * DVy2_Sh + 2 * DVy3_Sh + DVy4_Sh)
            VZ_Sh = VZ_Sh + dt / 6 * (DVz1_Sh + 2 * DVz2_Sh + 2 * DVz3_Sh + DVz4_Sh)

            plSystem.ReplaceSystem(X, Y, Z, VX, VY, VZ, X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, Alpha, F_dv)

            drPlanets = [planet.DrawedPlanet for planet in plSystem.planets]
            drTraces = [planet.DrawedTrace for planet in plSystem.planets]
            #self.SpWidget.canvas.axes.axis('scaled')
            self.SpWidget.canvas.axes.set(xlim=[-2 * Side+X_Sh, 2 * Side+X_Sh], ylim=[-Side+Y_Sh, Side+Y_Sh])

            return drPlanets + drTraces + [plSystem.spaceShip.DrawedSpaceShip] \
                   + [plSystem.spaceShip.DrawedTrace]

        fig = self.SpWidget.canvas.figure

        animation = FuncAnimation(fig, NewPoints, interval=dt * 1000, blit=True)

        self.SpWidget.canvas.draw()

app = QApplication([])
window = SpaceWidget()
window.show()
app.exec_()


