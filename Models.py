import math
import numpy as np
import random as rd
import sympy as sp
import pprint
import time
import scipy.io as io
import pickle


class PlanetSystem():
    def __init__(self, planets):
        self.planets = planets

    def add_new_planet(self, planet):
        self.planets.append(planet)

    def add_spaceship(self, spaceShip):
        self.spaceShip = spaceShip

    def replace_system(self, KSI, ETA, ZETA, VKSI, VETA, VZETA, KSI_Sh, ETA_Sh, ZETA_Sh, VKSI_Sh, VETA_Sh, VZETA_Sh, Phi_Sh=0, F_curr=0):
        for planet, ksi, eta, zeta, vksi, veta, vzeta in zip(self.planets, KSI, ETA, ZETA, VKSI, VETA, VZETA):
            planet.replace(ksi, eta, zeta, vksi, veta, vzeta)
            planet.re_draw()
        if (self.spaceShip):
            self.spaceShip.replace(KSI_Sh, ETA_Sh, ZETA_Sh, VKSI_Sh, VETA_Sh, VZETA_Sh, Phi_Sh, F_curr)
            self.spaceShip.re_draw()

    def replace_system_without_draw(self, KSI, ETA, ZETA, VKSI, VETA, VZETA, KSI_Sh, ETA_Sh, ZETA_Sh, VKSI_Sh, VETA_Sh, VZETA_Sh, Phi_Sh=0, F_curr=0):
        for planet, ksi, eta, zeta, vksi, veta, vzeta in zip(self.planets, KSI, ETA, ZETA, VKSI, VETA, VZETA):
            planet.replace(ksi, eta, zeta, vksi, veta, vzeta)
        if (self.spaceShip):
            self.spaceShip.replace(KSI_Sh, ETA_Sh, ZETA_Sh, VKSI_Sh, VETA_Sh, VZETA_Sh, Phi_Sh, F_curr)


    def draw(self, axes):
        for planet in self.planets:
            planet.draw(axes)
        if (self.spaceShip):
            self.spaceShip.draw(axes)

    def get_move_equations(self, obnul):
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
        # DVKSI = [
        #     sum([
        #         (planet.k* (ksi - cur_ksi)) / (sp.sqrt((ksi - cur_ksi) ** 2 + (eta - cur_eta) ** 2 + (zeta - cur_zeta) ** 2) ** 3)
        #         for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
        #         if (ksi != cur_ksi)
        #     ])
        #     for cur_ksi, cur_eta, cur_zeta, current_planet in zip(KSI, ETA, ZETA, self.planets)
        # ]

        # DVETA = [
        #     sum([
        #         (planet.k* (eta - cur_eta)) / (sp.sqrt((ksi - cur_ksi) ** 2 + (eta - cur_eta) ** 2 + (zeta - cur_zeta) ** 2) ** 3)
        #         for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
        #         if (ksi != cur_ksi)
        #     ])
        #     for cur_ksi, cur_eta, cur_zeta, current_planet in zip(KSI, ETA, ZETA, self.planets)
        # ]
        # DVZETA = [
        #     sum([
        #         (planet.k*(zeta - cur_zeta)) / (sp.sqrt((ksi - cur_ksi) ** 2 + (eta - cur_eta) ** 2 + (zeta - cur_zeta) ** 2) ** 3)
        #         for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
        #         if (ksi != cur_ksi)
        #     ])
        #     for cur_ksi, cur_eta, cur_zeta, current_planet in zip(KSI, ETA, ZETA, self.planets)
        # ]

        DVKSI = 0
        DVETA = 0
        DVZETA = 0


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

            Fx_dv_vs_Earth = F_dv * VKSI_Sh /(sp.sqrt(VKSI_Sh**2 + VETA_Sh**2))  # Сила x двигателя направленная против земли
            Fy_dv_vs_Earth = F_dv * VETA_Sh /(sp.sqrt(VKSI_Sh**2 + VETA_Sh**2))  # Сила y двигателя направленная против земли
            if(obnul):
                Fx_dv_vs_Earth = 0
                Fy_dv_vs_Earth = 0
            print(f'[xF] Fx_dv_vs_Earth', Fx_dv_vs_Earth)
            print(f'[xF] Fy_dv_vs_Earth', Fy_dv_vs_Earth)


            DVKSI_Sh = sum([
                (planet.k * (ksi - KSI_Sh)) / (sp.sqrt((ksi - KSI_Sh) ** 2 + (eta - ETA_Sh) ** 2 + (zeta - ZETA_Sh) ** 2) ** 3)
                for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
            ]) + Fx_dv_vs_Earth

            DVETA_Sh = sum([
                (planet.k * (eta - ETA_Sh)) / (sp.sqrt((ksi - KSI_Sh) ** 2 + (eta - ETA_Sh) ** 2 + (zeta - ZETA_Sh) ** 2) ** 3)
                for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
            ]) + Fy_dv_vs_Earth

            DVZETA_Sh = sum([
                (planet.k * (zeta - ZETA_Sh)) /  (sp.sqrt((ksi - KSI_Sh) ** 2 + (eta - ETA_Sh) ** 2 + (zeta - ZETA_Sh) ** 2) ** 3)
                for ksi, eta, zeta, planet in zip(KSI, ETA, ZETA, self.planets)
            ]) #+ F_dv / self.spaceShip.m * sp.sin(Beta)

        self.SpaceShipMoveEquations = sp.lambdify(
            [KSI_Sh, ETA_Sh, ZETA_Sh, VKSI_Sh, VETA_Sh, VZETA_Sh, KSI, ETA, ZETA, VKSI, VETA, VZETA, F_dv, Alpha,Beta],
            [DKSI_Sh, DETA_Sh, DZETA_Sh, DVKSI_Sh, DVETA_Sh, DVZETA_Sh])


    def get_state_vectors(self):
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

    def replace(self, ksi, eta, zeta, vksi, veta, vzeta):
        self.ksi = ksi
        self.eta = eta
        self.zeta = zeta
        self.Vksi = vksi
        self.Veta = veta
        self.Vzeta = vzeta

        self.TraceKSI = np.append(self.TraceKSI, ksi)
        self.TraceETA = np.append(self.TraceETA, eta)
        self.TraceZETA = np.append(self.TraceZETA, zeta)    

    def draw(self, axes):
        self.DrawedPlanet = axes.plot(self.ksi, self.eta, self.zeta, marker='o',markersize=self.R*50,color=self.color)[0]
        self.DrawedTrace = axes.plot(self.TraceKSI, self.TraceETA,self.TraceZETA, linestyle = ':',color=self.color)[0]

    def re_draw(self):
        self.DrawedPlanet.set_data_3d(self.ksi, self.eta,self.zeta)
        self.DrawedTrace.set_data_3d(self.TraceKSI, self.TraceETA,self.TraceZETA)

class SpaceShip():
    def __init__(self, ksi0, eta0, zeta0, Vksi0, Veta0,Vzeta0, m, R, color, F_max, K_stop_engine):
        self.ksi0 = ksi0
        self.eta0 = eta0
        self.zeta0 = zeta0
        self.Vksi0 = Vksi0
        self.Veta0 = Veta0
        self.Vzeta0 = Vzeta0
        self.m = m
        self.R = R
        self.color = color
        self.K_stop_engine = K_stop_engine

        self.ksi = ksi0
        self.eta = eta0
        self.zeta = zeta0
        self.Vksi = Vksi0
        self.Veta = Veta0
        self.Vzeta = Vzeta0

        self.phi = 0
        self.F_dv = F_max
        self.F_curr = 0

        self.SpaceShipX = self.ksi 
        self.SpaceShipY = self.eta
        self.SpaceShipZ = self.zeta

        self.SpaceShipFlameX = self.R * np.array([0, 1])
        self.SpaceShipFlameY = self.R * np.array([0, 1])

        self.TraceKSI = np.array([self.ksi])
        self.TraceETA = np.array([self.eta])
        self.TraceZETA = np.array([self.zeta])

    def replace(self, ksi, eta,zeta, vksi, veta, vzeta, phi, F_curr):
        self.ksi = ksi
        self.eta = eta
        self.zeta = zeta
        self.Vksi = vksi
        self.Veta = veta
        self.Vzeta = vzeta
        self.phi = phi
        self.F_curr = F_curr
        self.TraceKSI = np.append(self.TraceKSI, ksi)
        self.TraceETA = np.append(self.TraceETA, eta)
        self.TraceZETA = np.append(self.TraceZETA, zeta)


    def draw(self, axes):
        self.DrawedSpaceShip = axes.plot(self.ksi, self.eta, self.zeta,marker='o',markersize=1, color=self.color)[0]
        self.DrawedSpaceShipFlame = axes.plot(self.ksi + self.SpaceShipFlameX, self.eta + self.SpaceShipFlameY, self.zeta, color='yellow')[0]
        self.DrawedTrace = axes.plot(self.TraceKSI, self.TraceETA, self.TraceZETA,':')[0]

    def re_draw(self):
        #RotSpaceShipKSI, RotSpaceShipETA = rot_2D(self.SpaceShipKSI, self.SpaceShipETA, self.phi)
        self.DrawedSpaceShip.set_data_3d(self.ksi, self.eta, self.zeta)
        self.DrawedTrace.set_data_3d(self.TraceKSI, self.TraceETA,self.TraceZETA)

        RotSpaceShipFlameX, RotSpaceShipFlameY = rot_2D(self.ksi, self.eta, self.phi)
        self.DrawedSpaceShipFlame.set_data_3d(self.ksi + RotSpaceShipFlameX, self.eta + RotSpaceShipFlameY, self.zeta)



def rot_2D(KSI,ETA,phi):
    RotKSI = KSI*np.cos(phi) - ETA*np.sin(phi)
    RotETA = KSI*np.sin(phi) + ETA*np.cos(phi)
    return RotKSI, RotETA


def rot_3D(X, Y, Z, phi):
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

def draw_the_space(axes):
    axes.fill([-100, 100, 100, - 100], [-100, - 100, 100, 100], 'black')
    nstars = 5000
    xstars = 200 * np.random.random(nstars) - 100
    ystars = 200 * np.random.random(nstars) - 100
    brstars = 0.5 + np.random.random(nstars) / 2
    sizestars = np.random.random(nstars)
    for i in range(nstars):
        axes.plot(xstars[i], ystars[i], xstars[i], marker='o', markersize=sizestars[i], color=[brstars[i], brstars[i], brstars[i]])


