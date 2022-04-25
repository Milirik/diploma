import json
import numpy as np
import pprint
import pickle
import math
import os
import time
import random as rd
import sympy as sp
import scipy.io as io


from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
from matplotlib.animation import FuncAnimation

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *

from Models import PlanetSystem, Planet, SpaceShip
import FormOfSpaceObjects
import spacewidget


class SpaceWidget(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

class SpaceWidget(QMainWindow, FormOfSpaceObjects.Ui_MainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.koeffff = 0
        self.flag = False



        self.setupUi(self)
        self.setWindowTitle("Творение")

        self.fileData = []; # Данные из выбранного файла
        self.thisFile = ''; # Путь к выбранному файлу

        self.listEtudes.currentItemChanged.connect(self.chose_etude)
        self.createEtudeButton.clicked.connect(self.create_new_etude)
        self.deleteEtudeButton.clicked.connect(self.delete_etude)

        self.StartButton.clicked.connect(self.start_simulation)
        self.StartButton.clicked.connect(self.HereAreWeGo)          # ? так нужно...

        self.StopButton.clicked.connect(self.stop_simulation)

        self.listUniverseObjectsWidget.currentItemChanged.connect(self.chose_object)
        self.createObjectButton.clicked.connect(self.create_new_object)
        self.deleteObjectButton.clicked.connect(self.delete_object)
        self.saveDataButton.clicked.connect(self.save_changes_object)
        self.DrawSpaceshipTrajectory.clicked.connect(self.draw_spaceship_trajectory)


        self.get_all_etudes()
        self.typeObjectComboBox.addItems(['planet', 'starship'])

        self.addToolBar(NavigationToolbar(self.SpWidget.canvas, self))

        # self.AddPlanet.clicked.connect(self.AddPlanetFunction)
        # self.AddStar.clicked.connect(self.AddStarFunction)
        # self.AddShip.clicked.connect(self.AddShipFunction)
        # self.ShowSystem.clicked.connect(self.ShowSystemFunction)
        # self.deleteButton.clicked.connect(self.DeleteItemFunction)

    def get_all_etudes(self):
        self.listEtudes.clear()
        all_etudes = os.listdir('./etudes')
        self.listEtudes.addItems(all_etudes)

    def get_all_objects(self, file):
        fileData = []
        try:
            with open(file, "r") as read_file:
                fileData = json.load(read_file)
                self.fileData = fileData

            self.listUniverseObjectsWidget.clear()
            objects = [i['name']+' : '+i['type'] for i in fileData]
            self.listUniverseObjectsWidget.addItems(objects)
        except:
            print('[x] Error')

    def chose_etude(self):
        currentIndex = self.listEtudes.currentRow()
        if(currentIndex != -1):
            item = self.listEtudes.item(currentIndex).text()
            self.chosenEtudeLabel.setText(item)
            self.thisFile = os.path.abspath('.\\etudes\\' + item) 
            self.get_all_objects(self.thisFile)
            
    def create_new_etude(self):
        name_value = f'tmp{time.time()}' if self.nameEtudeField.text() == "" else self.nameEtudeField.text()
        with open(f"./etudes/{name_value}.json", "w") as write_file:
            json.dump([], write_file)
        self.get_all_etudes()

    def delete_etude(self):
        currentIndex = self.listEtudes.currentRow()
        item = '.\\etudes\\' + self.listEtudes.item(currentIndex).text()
        thisFile = os.path.abspath(item) 
        try:
            os.remove(thisFile)
            self.get_all_etudes()
            self.listUniverseObjectsWidget.clear()
        except:
            print(f'[x] Файл {thisFile} не удалось удалить')


    def start_simulation(self):
        pass

    def stop_simulation(self):
        self.animation.event_source.stop()


    def chose_object(self):
        
        currentIndex = self.listUniverseObjectsWidget.currentRow()

        if(currentIndex != -1):
            item = self.listUniverseObjectsWidget.item(currentIndex).text()
            self.chosenObjectLabel.setText(item)

            object_ = self.fileData[currentIndex]

            self.M_field.setText(str(object_["m"]))
            self.R_field.setText(str(object_["R"]))
            self.Color_field.setText(str(object_["color"]))
            self.P_field.setText(str(object_["P"]))

            self.x0_field.setText(str(object_["x"]))
            self.y0_field.setText(str(object_["y"]))
            self.z0_field.setText(str(object_["z"]))

            self.Vx0_field.setText(str(object_["Vx"]))
            self.Vy0_field.setText(str(object_["Vy"]))
            self.Vz0_field.setText(str(object_["Vz"]))

            self.F_dv_field.setText(str(object_["F_dv"]))
            self.Phi_field.setText(str(object_["Phi"]))
            self.Flamecolor_field.setText(str(object_["Fl_color"]))
  
    def create_new_object(self):
        name_value = f'tmp{time.time()}' if self.nameObjectField.text() == "" else self.nameObjectField.text()
        self.fileData.append({
            "name": name_value,
            "type": self.typeObjectComboBox.currentText(), 
            "m": 1.0,
            "R": 1.0,
            "color": "black",
            "P": 0.0,
            "x": 0.0,
            "y": 0.0,
            "z": 0.0,
            "Vx": 0.0,
            "Vy": 0.0,
            "Vz": 0.0,
            "F_dv": 0.0,
            "Phi": 0.0,
            "Fl_color": "yellow"
        })
        with open(self.thisFile, "w") as write_file:
            json.dump(self.fileData, write_file)

        self.get_all_objects(self.thisFile)

    def delete_object(self):
        currentIndex = self.listUniverseObjectsWidget.currentRow()
        self.fileData.pop(currentIndex)
        with open(self.thisFile, "w") as write_file:
            json.dump(self.fileData, write_file)

        self.get_all_objects(self.thisFile)

    def save_changes_object(self):
        currentIndex = self.listUniverseObjectsWidget.currentRow()
        self.fileData[currentIndex]['m'] = float(self.M_field.text())
        self.fileData[currentIndex]['R'] = float(self.R_field.text())
        self.fileData[currentIndex]['color'] = self.Color_field.text()
        self.fileData[currentIndex]['P'] = float(self.P_field.text())
        self.fileData[currentIndex]['x'] = float(self.x0_field.text())
        self.fileData[currentIndex]['y'] = float(self.y0_field.text())
        self.fileData[currentIndex]['z'] = float(self.z0_field.text())
        self.fileData[currentIndex]['Vx'] = float(self.Vx0_field.text())
        self.fileData[currentIndex]['Vy'] = float(self.Vy0_field.text())
        self.fileData[currentIndex]['Vz'] = float(self.Vz0_field.text())
        self.fileData[currentIndex]['F_dv'] = float(self.F_dv_field.text())
        self.fileData[currentIndex]['Phi'] = float(self.Phi_field.text())
        self.fileData[currentIndex]['Fl_color'] = self.Flamecolor_field.text()
        with open(self.thisFile, "w") as write_file:
            json.dump(self.fileData, write_file)

        self.get_all_objects(self.thisFile)



    def draw_spaceship_trajectory(self):
        self.progressBarDrawingSpTr.setValue(0) 
        self.DrawSpaceshipTrajectory.setEnabled(False)

        def NewPoints(i, traj):
            global t, dt, plSystem, ksi, eta, zeta, Vksi, Veta, Vzeta, Dksi, Deta, Dzeta, DVksi, DVeta, DVzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, Dksi_Sh, Deta_Sh, Dzeta_Sh, DVksi_Sh, DVeta_Sh, DVzeta_Sh, F_dv, Alpha, Beta
            t += 36000*dt

            #Методом Рунге - Кутты
            Dksi1, Deta1,Dzeta1, DVksi1, DVeta1,DVzeta1 = plSystem.SpaceBodyMoveEquations(ksi, eta, zeta, Vksi, Veta,Vzeta)
            Dksi1_Sh, Deta1_Sh,Dzeta1_Sh, DVksi1_Sh, DVeta1_Sh,DVzeta1_Sh = plSystem.SpaceShipMoveEquations(ksi_Sh, eta_Sh,zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, ksi, eta, zeta,Vksi, Veta, Vzeta, F_dv, Alpha,Beta)

            Dksi1 = np.array(Dksi1)
            Deta1 = np.array(Deta1)
            Dzeta1 = np.array(Dzeta1)
            DVksi1 = np.array(DVksi1)
            DVeta1 = np.array(DVeta1)
            DVzeta1 = np.array(DVzeta1)

            Dksi1_Sh = np.array(Dksi1_Sh)
            Deta1_Sh = np.array(Deta1_Sh)
            Dzeta1_Sh = np.array(Dzeta1_Sh)
            DVksi1_Sh = np.array(DVksi1_Sh)
            DVeta1_Sh = np.array(DVeta1_Sh)
            DVzeta1_Sh = np.array(DVzeta1_Sh)


            Dksi2, Deta2, Dzeta2, DVksi2, DVeta2, DVzeta2 = plSystem.SpaceBodyMoveEquations(ksi+Dksi1/2*dt, eta+Deta1/2*dt, zeta+Dzeta1/2*dt, Vksi+DVksi1/2*dt, Veta+DVeta1/2*dt,Vzeta+DVzeta1/2*dt)
            Dksi2_Sh, Deta2_Sh,Dzeta2_Sh, DVksi2_Sh, DVeta2_Sh,DVzeta2_Sh = plSystem.SpaceShipMoveEquations(
                ksi_Sh+Dksi1_Sh/2*dt, eta_Sh+Deta1_Sh/2*dt, zeta_Sh+Dzeta1_Sh/2*dt, Vksi_Sh+DVksi1_Sh/2*dt, Veta_Sh+DVeta1_Sh/2*dt,Vzeta_Sh+DVzeta1_Sh/2*dt,
                ksi+Dksi1/2*dt, eta+Deta1/2*dt,zeta+Dzeta1/2*dt, Vksi+DVksi1/2*dt, Veta+DVeta1/2*dt, Vzeta+DVzeta1/2*dt, F_dv, Alpha,Beta)

            Dksi2 = np.array(Dksi2)
            Deta2 = np.array(Deta2)
            Dzeta2 = np.array(Dzeta2)
            DVksi2 = np.array(DVksi2)
            DVeta2 = np.array(DVeta2)
            DVzeta2 = np.array(DVzeta2)
            Dksi2_Sh = np.array(Dksi2_Sh)
            Deta2_Sh = np.array(Deta2_Sh)
            Dzeta2_Sh = np.array(Dzeta2_Sh)
            DVksi2_Sh = np.array(DVksi2_Sh)
            DVeta2_Sh = np.array(DVeta2_Sh)
            DVzeta2_Sh = np.array(DVzeta2_Sh)

            Dksi3, Deta3, Dzeta3, DVksi3, DVeta3, DVzeta3 = plSystem.SpaceBodyMoveEquations(ksi+Dksi2/2*dt, eta+Deta2/2*dt, zeta+Dzeta2/2*dt, Vksi+DVksi2/2*dt, Veta+DVeta2/2*dt,Vzeta+DVzeta2/2*dt)
            Dksi3_Sh, Deta3_Sh, Dzeta3_Sh, DVksi3_Sh, DVeta3_Sh, DVzeta3_Sh = plSystem.SpaceShipMoveEquations(
                ksi_Sh + Dksi2_Sh / 2 * dt, eta_Sh + Deta2_Sh / 2 * dt, zeta_Sh + Dzeta2_Sh / 2 * dt,
                Vksi_Sh + DVksi2_Sh / 2 * dt, Veta_Sh + DVeta2_Sh / 2 * dt, Vzeta_Sh + DVzeta2_Sh / 2 * dt,
                ksi + Dksi2 / 2 * dt, eta + Deta2 / 2 * dt, zeta + Dzeta2 / 2 * dt, Vksi + DVksi2 / 2 * dt,
                Veta + DVeta2 / 2 * dt, Vzeta + DVzeta2 / 2 * dt, F_dv, Alpha, Beta)

            Dksi3 = np.array(Dksi3)
            Deta3 = np.array(Deta3)
            Dzeta3 = np.array(Dzeta3)
            DVksi3 = np.array(DVksi3)
            DVeta3 = np.array(DVeta3)
            DVzeta3 = np.array(DVzeta3)
            Dksi3_Sh = np.array(Dksi3_Sh)
            Deta3_Sh = np.array(Deta3_Sh)
            Dzeta3_Sh = np.array(Dzeta3_Sh)
            DVksi3_Sh = np.array(DVksi3_Sh)
            DVeta3_Sh = np.array(DVeta3_Sh)
            DVzeta3_Sh = np.array(DVzeta3_Sh)

            Dksi4, Deta4, Dzeta4, DVksi4, DVeta4, DVzeta4  = plSystem.SpaceBodyMoveEquations(ksi+Dksi3/2*dt, eta+Deta3/2*dt, zeta+Dzeta3/2*dt, Vksi+DVksi3/2*dt, Veta+DVeta3/2*dt,Vzeta+DVzeta3/2*dt)
            Dksi4_Sh, Deta4_Sh, Dzeta4_Sh, DVksi4_Sh, DVeta4_Sh, DVzeta4_Sh = plSystem.SpaceShipMoveEquations(
                ksi_Sh + Dksi3_Sh / 2 * dt, eta_Sh + Deta3_Sh / 2 * dt, zeta_Sh + Dzeta3_Sh / 2 * dt,
                Vksi_Sh + DVksi3_Sh / 2 * dt, Veta_Sh + DVeta3_Sh / 2 * dt, Vzeta_Sh + DVzeta3_Sh / 2 * dt,
                ksi + Dksi3 / 2 * dt, eta + Deta3 / 2 * dt, zeta + Dzeta3 / 2 * dt, Vksi + DVksi3 / 2 * dt,
                Veta + DVeta3 / 2 * dt, Vzeta + DVzeta3 / 2 * dt, F_dv, Alpha, Beta)

            Dksi4 = np.array(Dksi4)
            Deta4 = np.array(Deta4)
            Dzeta4 = np.array(Dzeta4)
            DVksi4 = np.array(DVksi4)
            DVeta4 = np.array(DVeta4)
            DVzeta4 = np.array(DVzeta4)
            Dksi4_Sh = np.array(Dksi4_Sh)
            Deta4_Sh = np.array(Deta4_Sh)
            Dzeta4_Sh = np.array(Dzeta4_Sh)
            DVksi4_Sh = np.array(DVksi4_Sh)
            DVeta4_Sh = np.array(DVeta4_Sh)
            DVzeta4_Sh = np.array(DVzeta4_Sh)

            ksi = ksi + dt/6 * (Dksi1 + 2*Dksi2 + 2*Dksi3 + Dksi4)
            eta = eta + dt/6 * (Deta1 + 2*Deta2 + 2*Deta3 + Deta4)
            zeta = zeta + dt / 6 * (Dzeta1 + 2 * Dzeta2 + 2 * Dzeta3 + Dzeta4)

            Vksi = Vksi + dt/6 * (DVksi1 + 2*DVksi2 + 2*DVksi3 + DVksi4)
            Veta = Veta + dt/6 * (DVeta1 + 2*DVeta2 + 2*DVeta3 + DVeta4)
            Vzeta = Vzeta + dt / 6 * (DVzeta1 + 2 * DVzeta2 + 2 * DVzeta3 + DVzeta4)

            ksi_Sh = ksi_Sh + dt / 6 * (Dksi1_Sh + 2 * Dksi2_Sh + 2 * Dksi3_Sh + Dksi4_Sh)
            eta_Sh = eta_Sh + dt / 6 * (Deta1_Sh + 2 * Deta2_Sh + 2 * Deta3_Sh + Deta4_Sh)
            zeta_Sh = zeta_Sh + dt / 6 * (Dzeta1_Sh + 2 * Dzeta2_Sh + 2 * Dzeta3_Sh + Dzeta4_Sh)

            Vksi_Sh = Vksi_Sh + dt / 6 * (DVksi1_Sh + 2 * DVksi2_Sh + 2 * DVksi3_Sh + DVksi4_Sh)
            Veta_Sh = Veta_Sh + dt / 6 * (DVeta1_Sh + 2 * DVeta2_Sh + 2 * DVeta3_Sh + DVeta4_Sh)
            Vzeta_Sh = Vzeta_Sh + dt / 6 * (DVzeta1_Sh + 2 * DVzeta2_Sh + 2 * DVzeta3_Sh + DVzeta4_Sh)

            if(self.koeffff > 3500 and not self.flag):
                plSystem.get_move_equations(True)
                self.flag = True
            else:
                self.koeffff+=1
            print(f'[x] ', plSystem.spaceShip.ksi, plSystem.spaceShip.eta, plSystem.spaceShip.zeta)


            plSystem.replace_system_without_draw(ksi, eta, zeta, Vksi, Veta,Vzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh)


            traj.append([plSystem.spaceShip.ksi, plSystem.spaceShip.ksi, plSystem.spaceShip.ksi])


            # return  [plSystem.spaceShip.DrawedTrace]

        global t, dt, plSystem, ksi, eta, zeta, Vksi, Veta, Vzeta, Dksi, Deta, Dzeta, DVksi, DVeta, DVzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, Dksi_Sh, Deta_Sh, Dzeta_Sh, DVksi_Sh, DVeta_Sh, DVzeta_Sh, F_dv, Alpha, Beta
        t = 0.0

        #     Параметры массы
        dt = float(self.TStep_field.text())

        # Было задано до этого
        # F_max = plSystem.spaceShip.F_dv
        # F_dv = self.F_Bar.value()*F_max/100
        # Alpha = self.Angle_Bar.value()/360*6.28+1.57

        F_dv = 0 #2500 # Сила двигателя
        Alpha = 0 #360/24*(t+dt) # Направленнность
        Beta = 0


        razm = 4.216424392e7 # Для обезразмеривания
        koff = 7.29e-5 # Для обезразмеривания


        plSystem = PlanetSystem([])
        for i in self.fileData:
            if(i['type'] == 'planet'):
                ksi_, eta_, zeta_ = [i["x"] / razm, i["y"] / razm, i["z"] / razm]
                print('[!]', ksi_, eta_, zeta_)
                V_ksi, V_eta, V_zeta = [i["Vx"]  / (koff * razm), i["Vy"]  / (koff * razm), i["Vz"]  / (koff * razm)]
                R =  i["R"] / razm
                M = i["m"]
                color = i["color"]

                ki = 0.9999999998 if i["name"] == 'Earth' else 0.01232376679 # Для обезразмеривания

                plSystem.add_new_planet(Planet(ksi_, eta_, zeta_, V_ksi, V_eta, V_zeta, ki, M, R, color))
            else:
                ksi_, eta_, zeta_ = [i["x"] / razm, i["y"] / razm, i["z"] / razm]
                V_ksi, V_eta, V_zeta = [i["Vx"]  / (koff * razm), i["Vy"]  / (koff * razm), i["Vz"]  / (koff * razm)]
                R =  6 * razm / razm
                M = i["m"]
                F_dv =  i["F_dv"]

                plSystem.add_spaceship(SpaceShip(ksi_, eta_, zeta_, V_ksi, V_eta, V_zeta, M, R, color, F_dv))

        # F_dv = 0 # Если убрать то будет норм двигатель работать для ракеты

        # ===================== Просчитываем траектории полета и получаем вектора ======================= #
        if((len(plSystem.planets) > 0 and hasattr(plSystem, "spaceShip")) or True): # Убрать TRUE
            plSystem.get_move_equations(self.koeffff)
            ksi,eta,zeta, Vksi, Veta,Vzeta = plSystem.get_state_vectors()
            ksi_Sh = plSystem.spaceShip.ksi
            eta_Sh = plSystem.spaceShip.eta
            zeta_Sh = plSystem.spaceShip.zeta
            Vksi_Sh = plSystem.spaceShip.Vksi
            Veta_Sh = plSystem.spaceShip.Veta
            Vzeta_Sh = plSystem.spaceShip.Vzeta 

        dt = 0.01
        cnt = 0
        max_cnt = 25000
        traj = []
        while cnt != max_cnt:
            NewPoints(dt, traj)
            cnt+=1
            self.progressBarDrawingSpTr.setValue(cnt*100/max_cnt) 

        self.DrawSpaceshipTrajectory.setEnabled(True)

        Side = float(self.K_field.text()) # Сторона графика. С помощью нее можно увеличить графики
        self.SpWidget.canvas.axes.clear()
        self.SpWidget.canvas.axes.set(xlim=[-2*Side, 2*Side], ylim=[-Side, Side], zlim=[-Side, Side])
        self.SpWidget.canvas.axes.set_title('Это космос')
        self.SpWidget.canvas.axes.set_xlabel('X')
        self.SpWidget.canvas.axes.set_ylabel('Y')
        self.SpWidget.canvas.axes.set_zlabel('Z')

        self.SpWidget.canvas.axes.plot(plSystem.spaceShip.TraceKSI, plSystem.spaceShip.TraceETA, plSystem.spaceShip.TraceZETA, ':')

        self.SpWidget.canvas.show()
        fig = self.SpWidget.canvas.figure
        self.SpWidget.canvas.draw()
   
    def HereAreWeGo(self):
        def NewPoints(i):
            global t, dt, plSystem, ksi, eta, zeta, Vksi, Veta, Vzeta, Dksi, Deta, Dzeta, DVksi, DVeta, DVzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, Dksi_Sh, Deta_Sh, Dzeta_Sh, DVksi_Sh, DVeta_Sh, DVzeta_Sh, F_dv, Alpha, Beta
            t += 36000*dt

            #Методом Рунге - Кутты
            Dksi1, Deta1,Dzeta1, DVksi1, DVeta1,DVzeta1 = plSystem.SpaceBodyMoveEquations(ksi, eta, zeta, Vksi, Veta,Vzeta)
            Dksi1_Sh, Deta1_Sh,Dzeta1_Sh, DVksi1_Sh, DVeta1_Sh,DVzeta1_Sh = plSystem.SpaceShipMoveEquations(ksi_Sh, eta_Sh,zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, ksi, eta, zeta,Vksi, Veta, Vzeta, F_dv, Alpha,Beta)

            #
            Dksi1 = np.array(Dksi1)
            Deta1 = np.array(Deta1)
            Dzeta1 = np.array(Dzeta1)
            DVksi1 = np.array(DVksi1)
            DVeta1 = np.array(DVeta1)
            DVzeta1 = np.array(DVzeta1)

            Dksi1_Sh = np.array(Dksi1_Sh)
            Deta1_Sh = np.array(Deta1_Sh)
            Dzeta1_Sh = np.array(Dzeta1_Sh)
            DVksi1_Sh = np.array(DVksi1_Sh)
            DVeta1_Sh = np.array(DVeta1_Sh)
            DVzeta1_Sh = np.array(DVzeta1_Sh)


            Dksi2, Deta2, Dzeta2, DVksi2, DVeta2, DVzeta2 = plSystem.SpaceBodyMoveEquations(ksi+Dksi1/2*dt, eta+Deta1/2*dt, zeta+Dzeta1/2*dt, Vksi+DVksi1/2*dt, Veta+DVeta1/2*dt,Vzeta+DVzeta1/2*dt)
            Dksi2_Sh, Deta2_Sh,Dzeta2_Sh, DVksi2_Sh, DVeta2_Sh,DVzeta2_Sh = plSystem.SpaceShipMoveEquations(
                ksi_Sh+Dksi1_Sh/2*dt, eta_Sh+Deta1_Sh/2*dt, zeta_Sh+Dzeta1_Sh/2*dt, Vksi_Sh+DVksi1_Sh/2*dt, Veta_Sh+DVeta1_Sh/2*dt,Vzeta_Sh+DVzeta1_Sh/2*dt,
                ksi+Dksi1/2*dt, eta+Deta1/2*dt,zeta+Dzeta1/2*dt, Vksi+DVksi1/2*dt, Veta+DVeta1/2*dt, Vzeta+DVzeta1/2*dt, F_dv, Alpha,Beta)

            Dksi2 = np.array(Dksi2)
            Deta2 = np.array(Deta2)
            Dzeta2 = np.array(Dzeta2)
            DVksi2 = np.array(DVksi2)
            DVeta2 = np.array(DVeta2)
            DVzeta2 = np.array(DVzeta2)
            Dksi2_Sh = np.array(Dksi2_Sh)
            Deta2_Sh = np.array(Deta2_Sh)
            Dzeta2_Sh = np.array(Dzeta2_Sh)
            DVksi2_Sh = np.array(DVksi2_Sh)
            DVeta2_Sh = np.array(DVeta2_Sh)
            DVzeta2_Sh = np.array(DVzeta2_Sh)

            Dksi3, Deta3, Dzeta3, DVksi3, DVeta3, DVzeta3 = plSystem.SpaceBodyMoveEquations(ksi+Dksi2/2*dt, eta+Deta2/2*dt, zeta+Dzeta2/2*dt, Vksi+DVksi2/2*dt, Veta+DVeta2/2*dt,Vzeta+DVzeta2/2*dt)
            Dksi3_Sh, Deta3_Sh, Dzeta3_Sh, DVksi3_Sh, DVeta3_Sh, DVzeta3_Sh = plSystem.SpaceShipMoveEquations(
                ksi_Sh + Dksi2_Sh / 2 * dt, eta_Sh + Deta2_Sh / 2 * dt, zeta_Sh + Dzeta2_Sh / 2 * dt,
                Vksi_Sh + DVksi2_Sh / 2 * dt, Veta_Sh + DVeta2_Sh / 2 * dt, Vzeta_Sh + DVzeta2_Sh / 2 * dt,
                ksi + Dksi2 / 2 * dt, eta + Deta2 / 2 * dt, zeta + Dzeta2 / 2 * dt, Vksi + DVksi2 / 2 * dt,
                Veta + DVeta2 / 2 * dt, Vzeta + DVzeta2 / 2 * dt, F_dv, Alpha, Beta)

            Dksi3 = np.array(Dksi3)
            Deta3 = np.array(Deta3)
            Dzeta3 = np.array(Dzeta3)
            DVksi3 = np.array(DVksi3)
            DVeta3 = np.array(DVeta3)
            DVzeta3 = np.array(DVzeta3)
            Dksi3_Sh = np.array(Dksi3_Sh)
            Deta3_Sh = np.array(Deta3_Sh)
            Dzeta3_Sh = np.array(Dzeta3_Sh)
            DVksi3_Sh = np.array(DVksi3_Sh)
            DVeta3_Sh = np.array(DVeta3_Sh)
            DVzeta3_Sh = np.array(DVzeta3_Sh)

            Dksi4, Deta4, Dzeta4, DVksi4, DVeta4, DVzeta4  = plSystem.SpaceBodyMoveEquations(ksi+Dksi3/2*dt, eta+Deta3/2*dt, zeta+Dzeta3/2*dt, Vksi+DVksi3/2*dt, Veta+DVeta3/2*dt,Vzeta+DVzeta3/2*dt)
            Dksi4_Sh, Deta4_Sh, Dzeta4_Sh, DVksi4_Sh, DVeta4_Sh, DVzeta4_Sh = plSystem.SpaceShipMoveEquations(
                ksi_Sh + Dksi3_Sh / 2 * dt, eta_Sh + Deta3_Sh / 2 * dt, zeta_Sh + Dzeta3_Sh / 2 * dt,
                Vksi_Sh + DVksi3_Sh / 2 * dt, Veta_Sh + DVeta3_Sh / 2 * dt, Vzeta_Sh + DVzeta3_Sh / 2 * dt,
                ksi + Dksi3 / 2 * dt, eta + Deta3 / 2 * dt, zeta + Dzeta3 / 2 * dt, Vksi + DVksi3 / 2 * dt,
                Veta + DVeta3 / 2 * dt, Vzeta + DVzeta3 / 2 * dt, F_dv, Alpha, Beta)

            Dksi4 = np.array(Dksi4)
            Deta4 = np.array(Deta4)
            Dzeta4 = np.array(Dzeta4)
            DVksi4 = np.array(DVksi4)
            DVeta4 = np.array(DVeta4)
            DVzeta4 = np.array(DVzeta4)
            Dksi4_Sh = np.array(Dksi4_Sh)
            Deta4_Sh = np.array(Deta4_Sh)
            Dzeta4_Sh = np.array(Dzeta4_Sh)
            DVksi4_Sh = np.array(DVksi4_Sh)
            DVeta4_Sh = np.array(DVeta4_Sh)
            DVzeta4_Sh = np.array(DVzeta4_Sh)

            ksi = ksi + dt/6 * (Dksi1 + 2*Dksi2 + 2*Dksi3 + Dksi4)
            eta = eta + dt/6 * (Deta1 + 2*Deta2 + 2*Deta3 + Deta4)
            zeta = zeta + dt / 6 * (Dzeta1 + 2 * Dzeta2 + 2 * Dzeta3 + Dzeta4)

            Vksi = Vksi + dt/6 * (DVksi1 + 2*DVksi2 + 2*DVksi3 + DVksi4)
            Veta = Veta + dt/6 * (DVeta1 + 2*DVeta2 + 2*DVeta3 + DVeta4)
            Vzeta = Vzeta + dt / 6 * (DVzeta1 + 2 * DVzeta2 + 2 * DVzeta3 + DVzeta4)

            ksi_Sh = ksi_Sh + dt / 6 * (Dksi1_Sh + 2 * Dksi2_Sh + 2 * Dksi3_Sh + Dksi4_Sh)
            eta_Sh = eta_Sh + dt / 6 * (Deta1_Sh + 2 * Deta2_Sh + 2 * Deta3_Sh + Deta4_Sh)
            zeta_Sh = zeta_Sh + dt / 6 * (Dzeta1_Sh + 2 * Dzeta2_Sh + 2 * Dzeta3_Sh + Dzeta4_Sh)

            Vksi_Sh = Vksi_Sh + dt / 6 * (DVksi1_Sh + 2 * DVksi2_Sh + 2 * DVksi3_Sh + DVksi4_Sh)
            Veta_Sh = Veta_Sh + dt / 6 * (DVeta1_Sh + 2 * DVeta2_Sh + 2 * DVeta3_Sh + DVeta4_Sh)
            Vzeta_Sh = Vzeta_Sh + dt / 6 * (DVzeta1_Sh + 2 * DVzeta2_Sh + 2 * DVzeta3_Sh + DVzeta4_Sh)

            if(self.koeffff > 3500 and not self.flag):
                plSystem.get_move_equations(True)
                self.flag = True
            else:
                self.koeffff+=1
            print(f'[x] ', ksi_Sh, eta_Sh, zeta_Sh)



            #print(ksi[1]-ksi[2], eta[1]-eta[2], zeta[1]-zeta[2])
            #print(sp.sqrt((ksi[1]-ksi[2])**2+(eta[1]-eta[2])**2+(zeta[1]-zeta[2])**2))
            plSystem.replace_system(ksi, eta, zeta, Vksi, Veta,Vzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh)

            #ctr_planet ='earth'
            #ax.set_xlim3d(ksi[1]- 165, ksi[1] + 165)
            #ax.set_ylim3d(eta[1]- 165, eta[1] + 165)
            #ax.set_zlim3d(zeta[1] - 165, zeta[1] + 165)

            drPlanets = [planet.DrawedPlanet for planet in plSystem.planets]
            drTraces = [planet.DrawedTrace for planet in plSystem.planets]

            return  [plSystem.spaceShip.DrawedSpaceShip]\
                   + drTraces+drPlanets+ [plSystem.spaceShip.DrawedTrace]

        global t, dt, plSystem, ksi, eta, zeta, Vksi, Veta, Vzeta, Dksi, Deta, Dzeta, DVksi, DVeta, DVzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, Dksi_Sh, Deta_Sh, Dzeta_Sh, DVksi_Sh, DVeta_Sh, DVzeta_Sh, F_dv, Alpha, Beta
        t = 0.0

        #     Параметры массы
        dt = float(self.TStep_field.text())

        # Было задано до этого
        # F_max = plSystem.spaceShip.F_dv
        # F_dv = self.F_Bar.value()*F_max/100
        # Alpha = self.Angle_Bar.value()/360*6.28+1.57

        F_dv = 0 #2500 # Сила двигателя
        Alpha = 0 #360/24*(t+dt) # Направленнность
        Beta = 0


        razm = 4.216424392e7 # Для обезразмеривания
        koff = 7.29e-5 # Для обезразмеривания



        plSystem = PlanetSystem([])
        for i in self.fileData:
            if(i['type'] == 'planet'):
                ksi_, eta_, zeta_ = [i["x"] / razm, i["y"] / razm, i["z"] / razm]
                print('[!]', ksi_, eta_, zeta_)

                V_ksi, V_eta, V_zeta = [i["Vx"]  / (koff * razm), i["Vy"]  / (koff * razm), i["Vz"]  / (koff * razm)]
                R =  i["R"] / razm
                M = i["m"]
                color = i["color"]

                ki = 0.9999999998 if i["name"] == 'Earth' else 0.01232376679 # Для обезразмеривания

                plSystem.add_new_planet(Planet(ksi_, eta_, zeta_, V_ksi, V_eta, V_zeta, ki, M, R, color))
            else:
                ksi_, eta_, zeta_ = [i["x"] / razm, i["y"] / razm, i["z"] / razm]
                V_ksi, V_eta, V_zeta = [i["Vx"]  / (koff * razm), i["Vy"]  / (koff * razm), i["Vz"]  / (koff * razm)]
                R =  6 * razm / razm
                M = i["m"]
                F_dv =  i["F_dv"]

                plSystem.add_spaceship(SpaceShip(ksi_, eta_, zeta_, V_ksi, V_eta, V_zeta, M, R, color, F_dv))

        # F_dv = 0 # Если убрать то будет норм двигатель работать для ракеты

        # ===================== Просчитываем траектории полета и получаем вектора ======================= #
        if((len(plSystem.planets) > 0 and hasattr(plSystem, "spaceShip")) or True): # Убрать TRUE
            plSystem.get_move_equations(self.koeffff)
            ksi,eta,zeta, Vksi, Veta,Vzeta = plSystem.get_state_vectors()
            ksi_Sh = plSystem.spaceShip.ksi
            eta_Sh = plSystem.spaceShip.eta
            zeta_Sh = plSystem.spaceShip.zeta
            Vksi_Sh = plSystem.spaceShip.Vksi
            Veta_Sh = plSystem.spaceShip.Veta
            Vzeta_Sh = plSystem.spaceShip.Vzeta


        # ====================================== Отрисовка графика ====================================== #
        Side = float(self.K_field.text()) # Сторона графика. С помощью нее можно увеличить графики
        # self.SpWidget.canvas.axes.clear()
        self.SpWidget.canvas.axes.set(xlim=[-2*Side, 2*Side], ylim=[-Side, Side], zlim=[-Side, Side])
        self.SpWidget.canvas.axes.set_title('Это космос')
        self.SpWidget.canvas.axes.set_xlabel('X')
        self.SpWidget.canvas.axes.set_ylabel('Y')
        self.SpWidget.canvas.axes.set_zlabel('Z')
        plSystem.draw(self.SpWidget.canvas.axes)
        self.SpWidget.canvas.show()
        fig = self.SpWidget.canvas.figure
        self.animation = FuncAnimation(fig, NewPoints, interval=dt * 1000, blit=True)
        self.SpWidget.canvas.draw()


if __name__ == "__main__":
    app = QApplication([])
    window = SpaceWidget()
    window.show()
    app.exec_()