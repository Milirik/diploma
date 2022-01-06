import json
import numpy as np
import pprint
import pickle
import FormOfSpaceObjects
import math
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *
import os
from UniverseSimulation import PlanetSystem, Planet, SpaceShip
import time

class SpaceWidget(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

class SpaceWidget(QMainWindow, FormOfSpaceObjects.Ui_MainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.setWindowTitle("Творение")

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
    

        # self.AddPlanet.clicked.connect(self.AddPlanetFunction)
        # self.AddStar.clicked.connect(self.AddStarFunction)
        # self.AddShip.clicked.connect(self.AddShipFunction)
        # self.ShowSystem.clicked.connect(self.ShowSystemFunction)
        # self.deleteButton.clicked.connect(self.DeleteItemFunction)

    def get_all_etudes(self):
        self.listEtudes.clear()
        all_etudes = os.listdir('./etudes')
        self.listEtudes.addItems(all_etudes)

    def chose_etude(self):
        currentIndex = self.listEtudes.currentRow()
        if(currentIndex != -1):
            item = self.listEtudes.item(currentIndex).text()
            self.chosenEtudeLabel.setText(item)
            thisFile = os.path.abspath('.\\etudes\\' + item) 

            fileData = []
            with open(thisFile, "r") as read_file:
                fileData = json.load(read_file)
                self.fileData = fileData

            self.listUniverseObjectsWidget.clear()
            objects = [i['name']+' : '+i['type'] for i in fileData]
            self.listUniverseObjectsWidget.addItems(objects)



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
        except:
            print(f'[x] Файл {thisFile} не удалось удалить')


    def start_simulation(self):
        pass

    def stop_simulation(self):
        pass


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
        pass

    def delete_object(self):
        pass

    def save_changes_object(self):
        currentIndex = self.listUniverseObjectsWidget.currentRow()
        item = self.listUniverseObjectsWidget.item(currentIndex).text()

        FileName = self.Fname_field.text() + '.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']

            print('!!!!!', item.split(' '))
            if('Космический' in item.split(' ')[1]):
                plSystem.spaceShip.x0 = float(self.x0_field.text())
                plSystem.spaceShip.y0 = float(self.y0_field.text())
                plSystem.spaceShip.z0 = float(self.z0_field.text())
                plSystem.spaceShip.Vx0 = float(self.Vx0_field.text())
                plSystem.spaceShip.Vy0 = float(self.Vy0_field.text())
                plSystem.spaceShip.Vz0 = float(self.Vz0_field.text())
                plSystem.spaceShip.m = float(self.M_field.text())
                plSystem.spaceShip.R = float(self.R_field.text())
                plSystem.spaceShip.color = [float(i) for i in self.Color_field.text().split(', ')]
                plSystem.spaceShip.flame_color = [float(i) for i in self.Flamecolor_field.text().split(', ')]

            else:
                plSystem.planets[currentIndex].x0 = float(self.x0_field.text())
                plSystem.planets[currentIndex].y0 = float(self.y0_field.text())
                plSystem.planets[currentIndex].z0 = float(self.z0_field.text())
                plSystem.planets[currentIndex].Vx0 = float(self.Vx0_field.text())
                plSystem.planets[currentIndex].Vy0 = float(self.Vy0_field.text())
                plSystem.planets[currentIndex].Vz0 = float(self.Vz0_field.text())
                plSystem.planets[currentIndex].m = float(self.M_field.text())
                plSystem.planets[currentIndex].R = float(self.R_field.text())
                plSystem.planets[currentIndex].color = [float(i) for i in self.Color_field.text().split(', ')]


            PS_dictionary = {'PS': plSystem}
            with open(FileName, 'wb') as PlanetSystemFile:
                pickle.dump(PS_dictionary, PlanetSystemFile)

            i = 1 
            self.listUniverseObjectsWidget.clear()
            objects = []
            for planet in plSystem.planets:
                str = f' \nПланета {i}: x0={planet.x0}, y0={planet.y0}, z0={planet.z0}, Vx0={planet.Vx0}, Vy0={planet.Vy0}, Vz0={planet.Vz0},  \nm={planet.m}, R={planet.R}, Цвет: {planet.color}\n'
                i += 1
                objects.append(str)
            if plSystem.spaceShip != []:
                str = f' \nКосмический корабль: \nx0={plSystem.spaceShip.x0}, y0={plSystem.spaceShip.y0}, z0={plSystem.spaceShip.z0},' \
                          f'Vx0={plSystem.spaceShip.Vx0}, Vy0={plSystem.spaceShip.Vy0}, Vz0={plSystem.spaceShip.Vz0}, \nm={plSystem.spaceShip.m}, ' \
                          f'R={plSystem.spaceShip.R}, Цвет: {plSystem.spaceShip.color}\n'
                objects.append(str)

            self.listUniverseObjectsWidget.addItems(objects)




   


    def DeleteItemFunction(self):
        currentIndex = self.listUniverseObjectsWidget.currentRow()
        item = self.listUniverseObjectsWidget.item(currentIndex).text()

        FileName = self.Fname_field.text() + '.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']

            if('Космический' in item.split(' ')[1]):
                plSystem.spaceShip = []
            elif('Планета' in item.split(' ')[1]):
                plSystem.planets.pop(currentIndex)
               
            PS_dictionary = {'PS': plSystem}
            with open(FileName, 'wb') as PlanetSystemFile:
                pickle.dump(PS_dictionary, PlanetSystemFile)

            self.ShowSystemFunction()
    
    def AddPlanetFunction(self):

        m = float(self.M_field.text())
        R = float(self.R_field.text())

        Color = [float(n) for n in self.Color_field.text().split(', ')]

        x = float(self.x0_field.text())
        y = float(self.y0_field.text())
        z = float(self.z0_field.text())
        Vx = float(self.Vx0_field.text())
        Vy = float(self.Vy0_field.text())
        Vz = float(self.Vz0_field.text())

        P = Planet(x, y, z, Vx, Vy, Vz, m, R, Color)

        FileName = self.Fname_field.text()+'.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']
            plSystem.AddNewPlanet(P)
        else:
            plSystem = PlanetSystem([P])
        PS_dictionary = {'PS': plSystem}
        with open(FileName, 'wb') as PlanetSystemFile:
            pickle.dump(PS_dictionary, PlanetSystemFile)
        self.ShowSystemFunction()

    def AddStarFunction(self):
        FileName = self.Fname_field.text()

    def AddShipFunction(self):

        m = float(self.M_field.text())
        R = float(self.R_field.text())

        Color = [float(n) for n in self.Color_field.text().split(', ')]

        x = float(self.x0_field.text())
        y = float(self.y0_field.text())
        z = float(self.z0_field.text())
        Vx = float(self.Vx0_field.text())
        Vy = float(self.Vy0_field.text())
        Vz = float(self.Vz0_field.text())

        Phi = float(self.Phi_field.text())
        F_dv = float(self.F_dv_field.text())
        Flame_color = [float(i) for i in self.Flamecolor_field.text().split(', ')]

        Sh = SpaceShip(x, y, z, Vx, Vy, Vz, m, R, Color, Flame_color, Phi, F_dv)

        FileName = self.Fname_field.text()+'.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']
            plSystem.AddSpaceShip(Sh)
        else:
            plSystem = PlanetSystem([], Sh)
        PS_dictionary = {'PS': plSystem}
        with open(FileName, 'wb') as PlanetSystemFile:
            pickle.dump(PS_dictionary, PlanetSystemFile)
        self.ShowSystemFunction()

    def ShowSystemFunction(self):
        FileName = self.Fname_field.text() + '.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']
            text = 'Планеты:'
            i = 1
            for planet in plSystem.planets:
                str = f' \nПланета {i}: x0={planet.x0}, y0={planet.y0}, z0={planet.z0}, Vx0={planet.Vx0}, Vy0={planet.Vy0}, Vz0={planet.Vz0}, \nm={planet.m}, R={planet.R}, Цвет: {planet.color}'
                i += 1
                text = text + str
            if plSystem.spaceShip != []:
                str = f' \nКосмический корабль: \nx0={plSystem.spaceShip.x0}, y0={plSystem.spaceShip.y0}, z0={plSystem.spaceShip.z0},' \
                      f'Vx0={plSystem.spaceShip.Vx0}, Vy0={plSystem.spaceShip.Vy0}, Vz0={plSystem.spaceShip.Vz0}, \nm={plSystem.spaceShip.m}, ' \
                      f'R={plSystem.spaceShip.R}, Цвет: {plSystem.spaceShip.color}\n'
                text = text + str
            self.UniverseBrowser.setText(text)
        else:
            self.UniverseBrowser.setText('Такого файла не существует')

        i = 1 
        self.listUniverseObjectsWidget.clear()

        objects = []
        for planet in plSystem.planets:
            str = f' \nПланета {i}: x0={planet.x0}, y0={planet.y0}, z0={planet.z0}, Vx0={planet.Vx0}, Vy0={planet.Vy0}, Vz0={planet.Vz0}, \nm={planet.m}, R={planet.R}, Цвет: {planet.color}'
            i += 1
            objects.append(str)
        if plSystem.spaceShip != []:
            str = f' \nКосмический корабль: \nx0={plSystem.spaceShip.x0}, y0={plSystem.spaceShip.y0}, z0={plSystem.spaceShip.z0},' \
                      f'Vx0={plSystem.spaceShip.Vx0}, Vy0={plSystem.spaceShip.Vy0}, Vz0={plSystem.spaceShip.Vz0}, \nm={plSystem.spaceShip.m}, ' \
                      f'R={plSystem.spaceShip.R}, Цвет: {plSystem.spaceShip.color}'
            objects.append(str)

        self.listUniverseObjectsWidget.addItems(objects)


if __name__ == "__main__":
    app = QApplication([])
    window = SpaceWidget()
    window.show()
    app.exec_()