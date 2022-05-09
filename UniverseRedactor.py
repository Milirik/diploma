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

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *

from SpaceSystemModelling import SpaceSystemModelling
import FormOfSpaceObjects
import spacewidget



class SpaceWidget(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

class SpaceWidget(QMainWindow, FormOfSpaceObjects.Ui_MainWindow, SpaceSystemModelling):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.setWindowTitle("Творение")

        self.fileData = []; # Данные из выбранного файла
        self.thisFile = ''; # Путь к выбранному файлу

        self.listEtudes.currentItemChanged.connect(self.chose_etude)
        self.createEtudeButton.clicked.connect(self.create_new_etude)
        self.deleteEtudeButton.clicked.connect(self.delete_etude)

        self.DrawSpaceshipTrajectory.clicked.connect(lambda: self.start_simulation(is_draw_only_trajectory=True))
        self.StartButton.clicked.connect(self.start_simulation)
        self.StopButton.clicked.connect(self.stop_simulation)
        self.ClearGraph.clicked.connect(self.clear_graph)

        self.listUniverseObjectsWidget.currentItemChanged.connect(self.chose_object)
        self.createObjectButton.clicked.connect(self.create_new_object)
        self.deleteObjectButton.clicked.connect(self.delete_object)
        self.saveDataButton.clicked.connect(self.save_changes_object)

        self.get_all_etudes()
        self.typeObjectComboBox.addItems(['planet', 'starship'])
        self.addToolBar(NavigationToolbar(self.SpWidget.canvas, self))


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

    def start_simulation(self, is_draw_only_trajectory):
        if(is_draw_only_trajectory): self.HereAreWeGo(is_draw_only_trajectory)
        else: self.HereAreWeGo()

    def stop_simulation(self):
        self.animation.event_source.stop()

    def clear_graph(self):
        self.SpWidget.canvas.axes.clear()
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
            self.K_stop_engine.setText(str(object_["K_stop_engine"]))
  
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
            "Fl_color": "yellow",
            "K_stop_engine": 0
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
        self.fileData[currentIndex]['K_stop_engine'] = int(self.K_stop_engine.text())

        with open(self.thisFile, "w") as write_file:
            json.dump(self.fileData, write_file)

        self.get_all_objects(self.thisFile)


if __name__ == "__main__":
    app = QApplication([])
    window = SpaceWidget()
    window.show()
    app.exec_()