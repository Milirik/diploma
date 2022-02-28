import math
import numpy as np
import random as rd
import sympy as sp
import pprint
import time
import scipy.io as io
import pickle


class PlanetSystem():
    def __init__(self, planets, spaceShip='NotEnoughGold'):
        self.planets = planets
        if spaceShip != 'NotEnoughGold':
            self.spaceShip = spaceShip

    def AddNewPlanet(self, planet):
        self.planets.append(planet)

    def AddSpaceShip(self, spaceShip):
        self.spaceShip = spaceShip

    def ReplaceSystem(self, KSI, ETA, ZETA, VKSI, VETA, VZETA, KSI_Sh, ETA_Sh, ZETA_Sh, VKSI_Sh, VETA_Sh, VZETA_Sh, Phi_Sh=0, F_curr=0):
        for planet, ksi, eta, zeta, vksi, veta, vzeta in zip(self.planets, KSI, ETA, ZETA, VKSI, VETA, VZETA):
            planet.Replace(ksi, eta, zeta, vksi, veta, vzeta)
            planet.ReDraw()
        if (self.spaceShip):
            self.spaceShip.Replace(KSI_Sh, ETA_Sh, ZETA_Sh, VKSI_Sh, VETA_Sh, VZETA_Sh, Phi_Sh, F_curr)
            self.spaceShip.ReDraw()

    # def ReplaceSystem(self, X, Y, Z, VX, VY, VZ, X_Sh=0, Y_Sh=0, Z_Sh=0, VX_Sh=0, VY_Sh=0, VZ_Sh=0, Phi_Sh=0, F_curr=0):
    #     for planet, x, y, z, vx, vy, vz in zip(self.planets, X, Y, Z, VX, VY, VZ):
    #         planet.Replace(x, y, z, vx, vy, vz)
    #         planet.ReDraw()
    #     if (self.spaceShip):
    #         self.spaceShip.Replace(X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, Phi_Sh, F_curr)
    #         self.spaceShip.ReDraw()

    def Draw(self, axes):
        for planet in self.planets:
            planet.Draw(axes)
        if (self.spaceShip):
            self.spaceShip.Draw(axes)

    def GetMoveEquations(self):
        n = len(self.planets)
        _strKSI = ''
        _strETA = ''
        _strZETA = ''
        _strVKSI = ''
        _strVETA = ''
        _strVZETA = ''
        for i in range(n):
            _strKSI += f'ksi{i}, '
            _strETA += f'eta{i}, '
            _strZETA += f'zeta{i}, '
            _strVKSI += f'Vksi{i}, '
            _strVETA += f'Veta{i}, '
            _strVZETA += f'Vzeta{i}, '

        G=6.674300000e-11

        KSI = sp.symbols(_strKSI)
        ETA = sp.symbols(_strETA)
        ZETA = sp.symbols(_strZETA)
        VKSI = sp.symbols(_strVKSI)
        VETA = sp.symbols(_strVETA)
        VZETA = sp.symbols(_strVZETA)

        DKSI= [Vksi for Vksi in VKSI]
        DETA = [Veta for Veta in VETA]
        DZETA = [Vzeta for Vzeta in VZETA]
        DVKSI = [
            sum([
                (planet.k* (ksi - cur_ksi)) / (sp.sqrt((ksi - cur_ksi) ** 2 + (eta - cur_eta) ** 2 + (zeta - cur_zeta) ** 2) ** 3)
                for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
                if (ksi != cur_ksi)
            ])
            for cur_ksi, cur_eta, cur_zeta, current_planet in zip(KSI, ETA, ZETA, self.planets)
        ]

        DVETA = [
            sum([
                (planet.k* (eta - cur_eta)) / (sp.sqrt((ksi - cur_ksi) ** 2 + (eta - cur_eta) ** 2 + (zeta - cur_zeta) ** 2) ** 3)
                for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
                if (ksi != cur_ksi)
            ])
            for cur_ksi, cur_eta, cur_zeta, current_planet in zip(KSI, ETA, ZETA, self.planets)
        ]
        DVZETA = [
            sum([
                (planet.k*(zeta - cur_zeta)) / (sp.sqrt((ksi - cur_ksi) ** 2 + (eta - cur_eta) ** 2 + (zeta - cur_zeta) ** 2) ** 3)
                for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
                if (ksi != cur_ksi)
            ])
            for cur_ksi, cur_eta, cur_zeta, current_planet in zip(KSI, ETA, ZETA, self.planets)
        ]
        self.SpaceBodyMoveEquations = sp.lambdify([KSI, ETA, ZETA, VKSI, VETA, VZETA], [DKSI, DETA, DZETA, DVKSI, DVETA, DVZETA])

        if (self.spaceShip):
            KSI_Sh = sp.symbols('ksi_Sh')
            ETA_Sh = sp.symbols('eta_Sh')
            ZETA_Sh = sp.symbols('zeta_Sh')
            VKSI_Sh = sp.symbols('Vksi_Sh')
            VETA_Sh = sp.symbols('Veta_Sh')
            VZETA_Sh = sp.symbols('Vzeta_Sh')

            F_dv = sp.symbols('f_dv')
            Alpha = sp.symbols('alpha')
            Beta = sp.symbols('beta')


            DKSI_Sh =VKSI_Sh
            DETA_Sh =VETA_Sh
            DZETA_Sh = VZETA_Sh

            DVKSI_Sh = sum([
                (planet.k* (ksi - KSI_Sh)) / (sp.sqrt((ksi - KSI_Sh) ** 2 + (eta - ETA_Sh) ** 2 + (zeta - ZETA_Sh) ** 2) ** 3)
                for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
            ]) + F_dv / self.spaceShip.m * sp.cos(Alpha)#*sp.cos(Beta)

            DVETA_Sh = sum([
                (planet.k* (eta - ETA_Sh)) / (sp.sqrt((ksi - KSI_Sh) ** 2 + (eta - ETA_Sh) ** 2 + (zeta - ZETA_Sh) ** 2) ** 3)
                for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
            ]) + F_dv / self.spaceShip.m * sp.sin(Alpha)#*sp.cos(Beta)
            DVZETA_Sh = sum([
                (planet.k* (zeta - ZETA_Sh)) /  (sp.sqrt((ksi - KSI_Sh) ** 2 + (eta - ETA_Sh) ** 2 + (zeta - ZETA_Sh) ** 2) ** 3)
                for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
            ]) #+ F_dv / self.spaceShip.m * sp.sin(Beta)

        self.SpaceShipMoveEquations = sp.lambdify(
            [KSI_Sh, ETA_Sh, ZETA_Sh, VKSI_Sh, VETA_Sh, VZETA_Sh, KSI, ETA, ZETA, VKSI, VETA, VZETA, F_dv, Alpha,Beta],
            [DKSI_Sh, DETA_Sh, DZETA_Sh, DVKSI_Sh, DVETA_Sh, DVZETA_Sh])

    # def GetMoveEquations(self):
    #     #    Земля                  Луна              Солнце    
    #     # R  6371км(3.6)            1738км(1)         696*10^3км(392.4)
    #     # m  5.9722*10^24кг(81.5)   7.35*10^22кг(1)   2*10^30кг(27139500)
    #     # 
    #     # Расстояния
    #     # от Земли до солнца = 149 600 000км (389.17)
    #     # от Земли до луны = 384 400км (1)
    #     #
    #     # Скорости
    #     # скорость вращения Земли вокруг Солнца = 29.765 км/с
    #     # скорость вращения Луны вокруг Земли 1.023 км/с
    #     #
    #     # t = Tau/omega
    #     # omega = 29.765 км/с

        

    #     n = len(self.planets)
    #     print(n, self.planets)
    #     _strX = ''
    #     _strY = ''
    #     _strZ = ''
    #     _strVx = ''
    #     _strVy = ''
    #     _strVz = ''
    #     for i in range(n):
    #         _strX += f'x{i}, '
    #         _strY += f'y{i}, '
    #         _strZ += f'z{i}, '
    #         _strVx += f'Vx{i}, '
    #         _strVy += f'Vy{i}, '
    #         _strVz += f'Vz{i}, '

    #     print(_strX)

    #     X = sp.symbols(_strX)
    #     Y = sp.symbols(_strY)
    #     Z = sp.symbols(_strZ)

    #     VX = sp.symbols(_strVx)
    #     VY = sp.symbols(_strVy)
    #     VZ = sp.symbols(_strVz)


    #     DX = [Vx for Vx in VX]
    #     DY = [Vy for Vy in VY]
    #     DZ = [Vz for Vz in VZ]

    #     DVX = [
    #         sum([
    #             (planet.m * (x - cur_x)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2 + (z - cur_z)**2) ** 3)
    #             for x, y, z, planet in zip(X, Y, Z, self.planets)
    #             if (x != cur_x)
    #         ])
    #         for cur_x, cur_y, cur_z, current_planet in zip(X, Y, Z, self.planets)
    #     ]

    #     DVY = [
    #         sum([
    #             (planet.m * (y - cur_y)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2 + (z - cur_z)**2) ** 3)
    #             for x, y, z, planet in zip(X, Y, Z, self.planets)
    #             if (x != cur_x)
    #         ])
    #         for cur_x, cur_y, cur_z, current_planet in zip(X, Y, Z, self.planets)
    #     ]

    #     DVZ = [
    #         sum([
    #             (planet.m * (z - cur_z)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2 + (z - cur_z)**2) ** 3)
    #             for x, y, z, planet in zip(X, Y, Z, self.planets)
    #             if (x != cur_x)
    #         ])
    #         for cur_x, cur_y, cur_z, current_planet in zip(X, Y, Z, self.planets)
    #     ]

    #     self.SpaceBodyMoveEquations = sp.lambdify([X, Y, Z, VX, VY, VZ], [DX, DY, DZ, DVX, DVY, DVZ])

    #     if (self.spaceShip):
    #         X_Sh = sp.symbols('x_Sh')
    #         Y_Sh = sp.symbols('y_Sh')
    #         Z_Sh = sp.symbols('z_Sh')

    #         VX_Sh = sp.symbols('Vx_Sh')
    #         VY_Sh = sp.symbols('Vy_Sh')
    #         VZ_Sh = sp.symbols('Vz_Sh')

    #         F_dv = sp.symbols('f_dv')
    #         Alpha = sp.symbols('alpha')
    #         Z_boost = sp.symbols('z_boost')


    #         DX_Sh = VX_Sh
    #         DY_Sh = VY_Sh
    #         DZ_Sh = VZ_Sh

    #         DVX_Sh = sum([
    #             (planet.m * (x - X_Sh)) / (sp.sqrt((x - X_Sh) ** 2 + (y - Y_Sh) ** 2 + (z - Z_Sh) ** 2) ** 3)
    #             for x, y, z, planet in zip(X, Y, Z, self.planets)
    #         ]) + F_dv / self.spaceShip.m * sp.cos(Z_boost)* sp.cos(Alpha)

    #         DVY_Sh = sum([
    #             (planet.m * (y - Y_Sh)) / (sp.sqrt((x - X_Sh) ** 2 + (y - Y_Sh) ** 2 + (z - Z_Sh) ** 2) ** 3)
    #             for x, y, z, planet in zip(X, Y, Z, self.planets)
    #         ]) + F_dv / self.spaceShip.m * sp.cos(Z_boost)* sp.sin(Alpha)


    #         DVZ_Sh = sum([
    #             (planet.m * (z - Z_Sh)) / (sp.sqrt((x - X_Sh) ** 2 + (y - Y_Sh) ** 2 + (z - Z_Sh) ** 2) ** 3)
    #             for x, y, z, planet in zip(X, Y, Z, self.planets)
    #         ]) + F_dv / self.spaceShip.m * sp.sin(Z_boost)


    #     self.SpaceShipMoveEquations = sp.lambdify([X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, X, Y, Z, VX, VY, VZ, F_dv, Alpha, Z_boost],
    #                                               [DX_Sh, DY_Sh, DZ_Sh, DVX_Sh, DVY_Sh, DVZ_Sh])

    def GetStateVectors(self):
        KSI = np.zeros(len(self.planets))
        ETA = np.zeros(len(self.planets))
        ZETA = np.zeros(len(self.planets))
        VKSI = np.zeros(len(self.planets))
        VETA = np.zeros(len(self.planets))
        VZETA = np.zeros(len(self.planets))
        for i in range(len(self.planets)):
            KSI[i] = self.planets[i].ksi
            ETA[i] = self.planets[i].eta
            ZETA[i] = self.planets[i].zeta
            VKSI[i] = self.planets[i].Vksi
            VETA[i] = self.planets[i].Veta
            VZETA[i] = self.planets[i].Vzeta

        return KSI, ETA, ZETA, VKSI, VETA, VZETA

    # def GetStateVectors(self):
    #     X = np.zeros(len(self.planets))
    #     Y = np.zeros(len(self.planets))
    #     Z = np.zeros(len(self.planets))
    #     VX = np.zeros(len(self.planets))
    #     VY = np.zeros(len(self.planets))
    #     VZ = np.zeros(len(self.planets))

    #     for i in range(len(self.planets)):
    #         print('[0]', self.planets[i].x, self.planets[i].y, self.planets[i].z, self.planets[i].Vx, self.planets[i].Vy, self.planets[i].Vz, type(self.planets[i].x))
    #         X[i] = self.planets[i].x
    #         Y[i] = self.planets[i].y
    #         Z[i] = self.planets[i].z

    #         VX[i] = self.planets[i].Vx
    #         VY[i] = self.planets[i].Vy
    #         VZ[i] = self.planets[i].Vz

    #     return X, Y, Z, VX, VY, VZ

class Planet():
    def __init__(self, ksi0, eta0, zeta0, Vksi0, Veta0, Vzeta0, k, m, R, color):
        self.ksi0 = ksi0
        self.eta0 = eta0
        self.zeta0 = zeta0
        self.Vksi0 = Vksi0
        self.Veta0 = Veta0
        self.Vzeta0 = Vzeta0
        self.k = k
        self.m = m
        self.R = R
        self.color = color

        self.ksi = ksi0
        self.eta = eta0
        self.zeta = zeta0
        self.Vksi = Vksi0
        self.Veta = Veta0
        self.Vzeta = Vzeta0

        phi = np.linspace(0, 6.28, 20)
        self.PlanetKSI = self.R * np.sin(phi)
        self.PlanetETA = self.R * np.cos(phi)
        self.PlanetZETA = self.R

        self.TraceKSI = np.array([self.ksi])
        self.TraceETA = np.array([self.eta])
        self.TraceZETA = np.array([self.zeta])
    # def __init__(self, x0, y0, z0, Vx0, Vy0, Vz0, m, R, color):
    #     self.x0 = x0
    #     self.y0 = y0
    #     self.z0 = z0
    #     self.Vx0 = Vx0
    #     self.Vy0 =  Vy0
    #     self.Vz0 = Vz0
    #     self.m = m
    #     self.R = R
    #     self.color = color

    #     self.x = x0
    #     self.y = y0
    #     self.z = z0
    #     self.Vx = Vx0
    #     self.Vy = Vy0
    #     self.Vz = Vz0

    #     phi = np.linspace(0, 6.28, 20)
    #     # self.PlanetX = self.R * np.sin(phi)
    #     # self.PlanetY = self.R * np.cos(phi)
    #     # self.PlanetZ = self.R * np.

    #     self.TraceX = np.array([self.x])
    #     self.TraceY = np.array([self.y])
    #     self.TraceZ = np.array([self.z])

    def Replace(self, ksi, eta, zeta, vksi, veta, vzeta):
        #print('Replaced planet')
        self.ksi = ksi
        self.eta = eta
        self.zeta = zeta
        self.Vksi = vksi
        self.Veta = veta
        self.Vzeta = vzeta

        self.TraceKSI = np.append(self.TraceKSI, ksi)
        self.TraceETA = np.append(self.TraceETA, eta)
        self.TraceZETA = np.append(self.TraceZETA, zeta)    
    # def Replace(self, x, y, z, vx, vy, vz):
    #     self.x = x
    #     self.y = y
    #     self.z = z
    #     self.Vx = vx
    #     self.Vy = vy
    #     self.Vz = vz

    #     self.TraceX = np.append(self.TraceX, x)
    #     self.TraceY = np.append(self.TraceY, y)
    #     self.TraceZ = np.append(self.TraceZ, z)
    def Draw(self, axes):
        self.DrawedPlanet = axes.plot(self.ksi, self.eta, self.zeta, marker='o',markersize=self.R*50,color=self.color)[0]
        self.DrawedTrace = axes.plot(self.TraceKSI, self.TraceETA,self.TraceZETA, linestyle = 'solid',color=self.color)[0]


    # def Draw(self, axes):
    #     # print('x y z', self.x, self.y, self.z)
    #     # print('Vx Vy Vz', self.Vx, self.Vy, self.Vz)
    #     # print('planet R', self.R)
    #     self.DrawedPlanet = axes.plot(self.x, self.y, self.z, color=self.color, marker='o', markersize=self.R)[0]
    #     self.DrawedTrace = axes.plot(self.TraceX, self.TraceY, self.TraceZ, ':')[0]

    def ReDraw(self):
        #print('Redrawedp planet')
        self.DrawedPlanet.set_data_3d(self.ksi, self.eta,self.zeta)
        self.DrawedTrace.set_data_3d(self.TraceKSI, self.TraceETA,self.TraceZETA)

    # def ReDraw(self):
    #     self.DrawedPlanet.set_data_3d(self.x, self.y, self.z)
    #     self.DrawedTrace.set_data_3d(self.TraceX, self.TraceY, self.TraceZ)

class SpaceShip():
    def __init__(self, ksi0, eta0, zeta0, Vksi0, Veta0,Vzeta0, m, R, color, F_max):
        self.ksi0 = ksi0
        self.eta0 = eta0
        self.zeta0 = zeta0
        self.Vksi0 = Vksi0
        self.Veta0 = Veta0
        self.Vzeta0 = Vzeta0
        self.m = m
        self.R = R
        self.color = color

        self.ksi = ksi0
        self.eta = eta0
        self.zeta = zeta0
        self.Vksi = Vksi0
        self.Veta = Veta0
        self.Vzeta = Vzeta0

        self.phi = 0
        self.F_dv = F_max
        self.F_curr = 0

        self.TraceKSI = np.array([self.ksi])
        self.TraceETA = np.array([self.eta])
        self.TraceZETA = np.array([self.zeta])
    # def __init__(self, x0, y0, z0, Vx0, Vy0, Vz0, m, R, color, FlameColor, Phi, F_max):
    #     self.x0 = x0
    #     self.y0 = y0
    #     self.z0 = z0
    #     self.Vx0 = Vx0
    #     self.Vy0 = Vy0
    #     self.Vz0 = Vz0
    #     self.m = m
    #     self.R = R
    #     self.color = color
    #     self.flame_color = FlameColor
    #     self.Phi = Phi

    #     self.x = x0
    #     self.y = y0
    #     self.z = z0
    #     self.Vx = Vx0
    #     self.Vy = Vy0
    #     self.Vz = Vz0

    #     self.phi = 0
    #     self.F_dv = F_max
    #     self.F_curr = 0

    #     self.SpaceShipX = self.x 
    #     self.SpaceShipY = self.y
    #     self.SpaceShipZ = self.z

    #     self.SpaceShipFlameX = self.R * np.array([0, -3.5,-3,-4,-3,-3.5, 0])
    #     self.SpaceShipFlameY = self.R * np.array([0.2, 0.15, 0.1, 0, -0.1, -0.15, -0.2])

    #     self.TraceX = np.array([self.x])
    #     self.TraceY = np.array([self.y])
    #     self.TraceZ = np.array([self.z])
    def Replace(self, ksi, eta,zeta, vksi, veta, vzeta, phi, F_curr):
        #print('Replaced sh')
        self.ksi = ksi
        self.eta = eta
        self.zeta = zeta
        self.Vksi = vksi
        self.Veta = veta
        self.Vzeta = vzeta
        self.phi = phi
        self.F_curr = F_curr
        #print(self.ksi,self.eta,self.zeta)
        self.TraceKSI = np.append(self.TraceKSI, ksi)
        self.TraceETA = np.append(self.TraceETA, eta)
        self.TraceZETA = np.append(self.TraceZETA, zeta)

    # def Replace(self, x, y, z, vx, vy, vz, phi, F_curr):
    #     self.x = x
    #     self.y = y
    #     self.z = z
    #     self.Vx = vx
    #     self.Vy = vy
    #     self.Vz = vz
    #     self.phi = phi
    #     self.F_curr = F_curr

    #     self.TraceX = np.append(self.TraceX, x)
    #     self.TraceY = np.append(self.TraceY, y)
    #     self.TraceZ = np.append(self.TraceZ, z)
    def Draw(self, axes):
        self.DrawedSpaceShip = axes.plot(self.ksi, self.eta, self.zeta,marker='o',markersize=1, color=self.color)[0]
        self.DrawedTrace = axes.plot(self.TraceKSI, self.TraceETA, self.TraceZETA,':')[0]

    # def Draw(self, axes):
    #     print('spaceShip R', self.R)

    #     self.DrawedSpaceShip = axes.plot(self.x, self.y, self.z, color='black', marker='o', markersize=self.R)[0]
    #     # self.DrawedSpaceShipFlame = axes.plot(self.x + self.SpaceShipFlameX, self.y + self.SpaceShipFlameY, self.z, color=self.flame_color)[0]
    #     self.DrawedTrace = axes.plot(self.TraceX, self.TraceY, self.TraceZ, ':')[0]

    def ReDraw(self):
        #print('Redrawed sh')
        #RotSpaceShipKSI, RotSpaceShipETA = Rot2D(self.SpaceShipKSI, self.SpaceShipETA, self.phi)
        self.DrawedSpaceShip.set_data_3d(self.ksi, self.eta, self.zeta)
        self.DrawedTrace.set_data_3d(self.TraceKSI, self.TraceETA,self.TraceZETA)
    # def ReDraw(self):
    #     # RotSpaceShipX, RotSpaceShipY, RotSpaceShipZ = Rot3D(self.x, self.y, self.z, self.phi)
    #     self.DrawedSpaceShip.set_data_3d(self.x, self.y, self.z)

    #     # PFlameX = self.F_curr/self.F_dv*self.SpaceShipFlameX + self.SpaceShipX[int(len(self.SpaceShipX)/2)]
    #     # RotSpaceShipFlameX, RotSpaceShipFlameY, RotSpaceShipFlameZ= Rot3D(PFlameX, self.SpaceShipFlameY, self.phi)
    #     # self.DrawedSpaceShipFlame.set_data_3d(self.x + RotSpaceShipFlameX, self.y + RotSpaceShipFlameY, self.z + RotSpaceShipFlameZ)
    #     self.DrawedTrace.set_data_3d(self.TraceX, self.TraceY, self.TraceZ)

def Rot2D(KSI,ETA,phi): # не было у меня
    RotKSI = KSI*np.cos(phi) - ETA*np.sin(phi)
    RotETA = KSI*np.sin(phi) + ETA*np.cos(phi)
    return RotKSI, RotETA


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


