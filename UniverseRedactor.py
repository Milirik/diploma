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
        self.setupUi(self)
        self.setWindowTitle("Творение")

        self.fileData = []; # Данные из выбранного файла
        self.thisFile = ''; # Путь к выбранному файлу

        self.listEtudes.currentItemChanged.connect(self.chose_etude)
        self.createEtudeButton.clicked.connect(self.create_new_etude)
        self.deleteEtudeButton.clicked.connect(self.delete_etude)

        self.StartButton.clicked.connect(self.start_simulation)
        self.StopButton.clicked.connect(self.stop_simulation)

        self.listUniverseObjectsWidget.currentItemChanged.connect(self.chose_object)
        self.createObjectButton.clicked.connect(self.create_new_object)
        self.deleteObjectButton.clicked.connect(self.delete_object)
        self.saveDataButton.clicked.connect(self.save_changes_object)

        self.get_all_etudes()
        self.typeObjectComboBox.addItems(['planet', 'starship'])

        self.StartButton.clicked.connect(self.HereAreWeGo)
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

   
    def HereAreWeGo(self):
        def NewPoints(i):
            global t, dt, plSystem, X, Y, Z, VX, VY, VZ, Dx, Dy, Dz, DVx, DVy, DVz, X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, Dx_Sh, Dy_Sh, Dz_Sh, DVx_Sh, DVy_Sh, DVz_Sh, F_max, F_dv, Alpha
            t += dt
            F_dv = self.F_Bar.value()*F_max/100
            Alpha = -self.Angle_Bar.value()/360*6.28+1.57
            Z_boost = self.Z_turbo.value()/360*6.28

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
        global t, dt, plSystem, X, Y, Z, VX, VY, VZ, Dx, Dy, Dz, DVx, DVy, DVz, X_Sh, Y_Sh, Z_Sh, VX_Sh, VY_Sh, VZ_Sh, Dx_Sh, Dy_Sh, Dz_Sh, DVx_Sh, DVy_Sh, DVZ_Sh, F_max, F_dv, Alpha
        
        #     Параметры массы
        dt = float(self.TStep_field.text())
        K = float(self.K_field.text())

        plSystem = PlanetSystem([])

        for i in self.fileData:
            if(i['type'] == 'planet'):
                plSystem.AddNewPlanet(Planet(i["x"], i["y"], i["z"], i["Vx"], i["Vy"], i["Vz"], i["m"], i["R"], i["color"]))
            else:
                plSystem.AddSpaceShip(SpaceShip(i["x"], i["y"], i["z"], i["Vx"], i["Vy"], i["Vz"], i["m"], i["R"], i["color"], i["Fl_color"], i["Phi"], i["F_dv"]))

        print('[u]', len(plSystem.planets), hasattr(plSystem, "spaceShip"))
        if(len(plSystem.planets) > 0 and hasattr(plSystem, "spaceShip")):

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




            fig = self.SpWidget.canvas.figure

            self.animation = FuncAnimation(fig, NewPoints, interval=dt * 1000, blit=True)

            self.SpWidget.canvas.draw()

        

            

if __name__ == "__main__":
    app = QApplication([])
    window = SpaceWidget()
    window.show()
    app.exec_()