# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'FormOfSpaceObjects.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1400, 1000)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.listUniverseObjectsWidget = QtWidgets.QListWidget(self.centralwidget)
        self.listUniverseObjectsWidget.setGeometry(QtCore.QRect(940, 79, 161, 551))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.listUniverseObjectsWidget.setFont(font)
        self.listUniverseObjectsWidget.setObjectName("listUniverseObjectsWidget")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(60, 10, 80, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(22)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.label_23 = QtWidgets.QLabel(self.centralwidget)
        self.label_23.setGeometry(QtCore.QRect(1260, 350, 21, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_23.setFont(font)
        self.label_23.setObjectName("label_23")
        self.Color = QtWidgets.QLabel(self.centralwidget)
        self.Color.setGeometry(QtCore.QRect(1130, 150, 51, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.Color.setFont(font)
        self.Color.setObjectName("Color")
        self.saveDataButton = QtWidgets.QPushButton(self.centralwidget)
        self.saveDataButton.setGeometry(QtCore.QRect(1230, 570, 151, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.saveDataButton.setFont(font)
        self.saveDataButton.setStyleSheet("background-color: rgb(85, 170, 255);")
        self.saveDataButton.setObjectName("saveDataButton")
        self.z0_field = QtWidgets.QLineEdit(self.centralwidget)
        self.z0_field.setGeometry(QtCore.QRect(1160, 350, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.z0_field.setFont(font)
        self.z0_field.setObjectName("z0_field")
        self.R_field = QtWidgets.QLineEdit(self.centralwidget)
        self.R_field.setGeometry(QtCore.QRect(1290, 110, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.R_field.setFont(font)
        self.R_field.setObjectName("R_field")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(1270, 110, 21, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.x0_field = QtWidgets.QLineEdit(self.centralwidget)
        self.x0_field.setGeometry(QtCore.QRect(1160, 270, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.x0_field.setFont(font)
        self.x0_field.setObjectName("x0_field")
        self.label_22 = QtWidgets.QLabel(self.centralwidget)
        self.label_22.setGeometry(QtCore.QRect(1190, 440, 31, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_22.setFont(font)
        self.label_22.setObjectName("label_22")
        self.Vy0_field = QtWidgets.QLineEdit(self.centralwidget)
        self.Vy0_field.setGeometry(QtCore.QRect(1290, 310, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.Vy0_field.setFont(font)
        self.Vy0_field.setObjectName("Vy0_field")
        self.F_dv_field = QtWidgets.QLineEdit(self.centralwidget)
        self.F_dv_field.setGeometry(QtCore.QRect(1240, 400, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.F_dv_field.setFont(font)
        self.F_dv_field.setObjectName("F_dv_field")
        self.y0_field = QtWidgets.QLineEdit(self.centralwidget)
        self.y0_field.setGeometry(QtCore.QRect(1160, 310, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.y0_field.setFont(font)
        self.y0_field.setObjectName("y0_field")
        self.Phi_field = QtWidgets.QLineEdit(self.centralwidget)
        self.Phi_field.setGeometry(QtCore.QRect(1240, 440, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.Phi_field.setFont(font)
        self.Phi_field.setObjectName("Phi_field")
        self.Color_field = QtWidgets.QLineEdit(self.centralwidget)
        self.Color_field.setGeometry(QtCore.QRect(1170, 150, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.Color_field.setFont(font)
        self.Color_field.setObjectName("Color_field")
        self.Color_2 = QtWidgets.QLabel(self.centralwidget)
        self.Color_2.setGeometry(QtCore.QRect(1270, 150, 31, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.Color_2.setFont(font)
        self.Color_2.setObjectName("Color_2")
        self.label_17 = QtWidgets.QLabel(self.centralwidget)
        self.label_17.setGeometry(QtCore.QRect(1130, 310, 21, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_17.setFont(font)
        self.label_17.setObjectName("label_17")
        self.P_field = QtWidgets.QLineEdit(self.centralwidget)
        self.P_field.setGeometry(QtCore.QRect(1290, 150, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.P_field.setFont(font)
        self.P_field.setObjectName("P_field")
        self.label_20 = QtWidgets.QLabel(self.centralwidget)
        self.label_20.setGeometry(QtCore.QRect(1260, 270, 21, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_20.setFont(font)
        self.label_20.setObjectName("label_20")
        self.label_18 = QtWidgets.QLabel(self.centralwidget)
        self.label_18.setGeometry(QtCore.QRect(1190, 400, 31, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_18.setFont(font)
        self.label_18.setObjectName("label_18")
        self.Flamecolor_field = QtWidgets.QLineEdit(self.centralwidget)
        self.Flamecolor_field.setGeometry(QtCore.QRect(1240, 480, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.Flamecolor_field.setFont(font)
        self.Flamecolor_field.setObjectName("Flamecolor_field")
        self.Vx0_field = QtWidgets.QLineEdit(self.centralwidget)
        self.Vx0_field.setGeometry(QtCore.QRect(1290, 270, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.Vx0_field.setFont(font)
        self.Vx0_field.setObjectName("Vx0_field")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(1130, 110, 21, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.label_21 = QtWidgets.QLabel(self.centralwidget)
        self.label_21.setGeometry(QtCore.QRect(1260, 310, 21, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_21.setFont(font)
        self.label_21.setObjectName("label_21")
        self.label_16 = QtWidgets.QLabel(self.centralwidget)
        self.label_16.setGeometry(QtCore.QRect(1130, 270, 21, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_16.setFont(font)
        self.label_16.setObjectName("label_16")
        self.Color_3 = QtWidgets.QLabel(self.centralwidget)
        self.Color_3.setGeometry(QtCore.QRect(1190, 480, 51, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.Color_3.setFont(font)
        self.Color_3.setObjectName("Color_3")
        self.Vz0_field = QtWidgets.QLineEdit(self.centralwidget)
        self.Vz0_field.setGeometry(QtCore.QRect(1290, 350, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.Vz0_field.setFont(font)
        self.Vz0_field.setObjectName("Vz0_field")
        self.label_19 = QtWidgets.QLabel(self.centralwidget)
        self.label_19.setGeometry(QtCore.QRect(1130, 350, 21, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_19.setFont(font)
        self.label_19.setObjectName("label_19")
        self.M_field = QtWidgets.QLineEdit(self.centralwidget)
        self.M_field.setGeometry(QtCore.QRect(1170, 110, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.M_field.setFont(font)
        self.M_field.setObjectName("M_field")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(960, 10, 111, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(22)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(1160, 10, 231, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(22)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(500, 10, 141, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(22)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.SpWidget = SpaceWidget(self.centralwidget)
        self.SpWidget.setGeometry(QtCore.QRect(230, 80, 671, 551))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.SpWidget.setFont(font)
        self.SpWidget.setObjectName("SpWidget")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(320, 750, 141, 20))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.TStep_field = QtWidgets.QLineEdit(self.centralwidget)
        self.TStep_field.setGeometry(QtCore.QRect(430, 750, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.TStep_field.setFont(font)
        self.TStep_field.setObjectName("TStep_field")
        self.StopButton = QtWidgets.QPushButton(self.centralwidget)
        self.StopButton.setGeometry(QtCore.QRect(540, 860, 121, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.StopButton.setFont(font)
        self.StopButton.setStyleSheet("background-color: rgb(255, 255, 127);")
        self.StopButton.setObjectName("StopButton")
        self.K_field = QtWidgets.QLineEdit(self.centralwidget)
        self.K_field.setGeometry(QtCore.QRect(430, 790, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.K_field.setFont(font)
        self.K_field.setObjectName("K_field")
        self.label_11 = QtWidgets.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(320, 790, 111, 20))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.StartButton = QtWidgets.QPushButton(self.centralwidget)
        self.StartButton.setGeometry(QtCore.QRect(400, 860, 121, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.StartButton.setFont(font)
        self.StartButton.setStyleSheet("background-color: rgb(85, 255, 127);")
        self.StartButton.setObjectName("StartButton")
        self.label_24 = QtWidgets.QLabel(self.centralwidget)
        self.label_24.setGeometry(QtCore.QRect(320, 710, 181, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setUnderline(True)
        self.label_24.setFont(font)
        self.label_24.setObjectName("label_24")
        self.label_25 = QtWidgets.QLabel(self.centralwidget)
        self.label_25.setGeometry(QtCore.QRect(30, 640, 151, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setUnderline(True)
        self.label_25.setFont(font)
        self.label_25.setObjectName("label_25")
        self.deleteLabel = QtWidgets.QLabel(self.centralwidget)
        self.deleteLabel.setGeometry(QtCore.QRect(30, 760, 141, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setUnderline(True)
        self.deleteLabel.setFont(font)
        self.deleteLabel.setObjectName("deleteLabel")
        self.listEtudes = QtWidgets.QListWidget(self.centralwidget)
        self.listEtudes.setGeometry(QtCore.QRect(30, 78, 150, 551))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.listEtudes.setFont(font)
        self.listEtudes.setObjectName("listEtudes")
        self.nameEtudeField = QtWidgets.QLineEdit(self.centralwidget)
        self.nameEtudeField.setGeometry(QtCore.QRect(90, 680, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.nameEtudeField.setFont(font)
        self.nameEtudeField.setText("")
        self.nameEtudeField.setObjectName("nameEtudeField")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(30, 680, 61, 20))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.createEtudeButton = QtWidgets.QPushButton(self.centralwidget)
        self.createEtudeButton.setGeometry(QtCore.QRect(90, 720, 91, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.createEtudeButton.setFont(font)
        self.createEtudeButton.setStyleSheet("background-color: rgb(85, 85, 127);\n"
"color: rgb(255, 255, 255);")
        self.createEtudeButton.setObjectName("createEtudeButton")
        self.deleteEtudeButton = QtWidgets.QPushButton(self.centralwidget)
        self.deleteEtudeButton.setGeometry(QtCore.QRect(30, 800, 151, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.deleteEtudeButton.setFont(font)
        self.deleteEtudeButton.setStyleSheet("background-color:#ff5e42;")
        self.deleteEtudeButton.setAutoDefault(False)
        self.deleteEtudeButton.setObjectName("deleteEtudeButton")
        self.label_28 = QtWidgets.QLabel(self.centralwidget)
        self.label_28.setGeometry(QtCore.QRect(940, 640, 161, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setUnderline(True)
        self.label_28.setFont(font)
        self.label_28.setObjectName("label_28")
        self.label_29 = QtWidgets.QLabel(self.centralwidget)
        self.label_29.setGeometry(QtCore.QRect(940, 790, 161, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setUnderline(True)
        self.label_29.setFont(font)
        self.label_29.setObjectName("label_29")
        self.deleteObjectButton = QtWidgets.QPushButton(self.centralwidget)
        self.deleteObjectButton.setGeometry(QtCore.QRect(940, 830, 161, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.deleteObjectButton.setFont(font)
        self.deleteObjectButton.setAutoFillBackground(False)
        self.deleteObjectButton.setStyleSheet("background-color: rgb(255, 94, 66);")
        self.deleteObjectButton.setObjectName("deleteObjectButton")
        self.label_30 = QtWidgets.QLabel(self.centralwidget)
        self.label_30.setGeometry(QtCore.QRect(1130, 70, 101, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setUnderline(True)
        self.label_30.setFont(font)
        self.label_30.setObjectName("label_30")
        self.label_31 = QtWidgets.QLabel(self.centralwidget)
        self.label_31.setGeometry(QtCore.QRect(1130, 220, 181, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setUnderline(True)
        self.label_31.setFont(font)
        self.label_31.setObjectName("label_31")
        self.createObjectButton = QtWidgets.QPushButton(self.centralwidget)
        self.createObjectButton.setGeometry(QtCore.QRect(1010, 750, 91, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.createObjectButton.setFont(font)
        self.createObjectButton.setStyleSheet("background-color: rgb(85, 85, 127);\n"
"color: rgb(255, 255, 255);")
        self.createObjectButton.setObjectName("createObjectButton")
        self.typeObjectComboBox = QtWidgets.QComboBox(self.centralwidget)
        self.typeObjectComboBox.setGeometry(QtCore.QRect(1010, 680, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        self.typeObjectComboBox.setFont(font)
        self.typeObjectComboBox.setObjectName("typeObjectComboBox")
        self.label_15 = QtWidgets.QLabel(self.centralwidget)
        self.label_15.setGeometry(QtCore.QRect(940, 680, 61, 20))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_15.setFont(font)
        self.label_15.setObjectName("label_15")
        self.nameObjectField = QtWidgets.QLineEdit(self.centralwidget)
        self.nameObjectField.setGeometry(QtCore.QRect(1010, 710, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.nameObjectField.setFont(font)
        self.nameObjectField.setText("")
        self.nameObjectField.setObjectName("nameObjectField")
        self.label_32 = QtWidgets.QLabel(self.centralwidget)
        self.label_32.setGeometry(QtCore.QRect(940, 710, 61, 20))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_32.setFont(font)
        self.label_32.setObjectName("label_32")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(30, 50, 47, 13))
        self.label_2.setObjectName("label_2")
        self.chosenEtudeLabel = QtWidgets.QLabel(self.centralwidget)
        self.chosenEtudeLabel.setGeometry(QtCore.QRect(80, 50, 101, 16))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.chosenEtudeLabel.setFont(font)
        self.chosenEtudeLabel.setText("")
        self.chosenEtudeLabel.setObjectName("chosenEtudeLabel")
        self.label_34 = QtWidgets.QLabel(self.centralwidget)
        self.label_34.setGeometry(QtCore.QRect(940, 50, 47, 13))
        self.label_34.setObjectName("label_34")
        self.chosenObjectLabel = QtWidgets.QLabel(self.centralwidget)
        self.chosenObjectLabel.setGeometry(QtCore.QRect(990, 50, 111, 16))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.chosenObjectLabel.setFont(font)
        self.chosenObjectLabel.setText("")
        self.chosenObjectLabel.setObjectName("chosenObjectLabel")
        self.DrawSpaceshipTrajectory = QtWidgets.QPushButton(self.centralwidget)
        self.DrawSpaceshipTrajectory.setEnabled(True)
        self.DrawSpaceshipTrajectory.setGeometry(QtCore.QRect(450, 650, 211, 23))
        self.DrawSpaceshipTrajectory.setObjectName("DrawSpaceshipTrajectory")
        self.progressBarDrawingSpTr = QtWidgets.QProgressBar(self.centralwidget)
        self.progressBarDrawingSpTr.setGeometry(QtCore.QRect(670, 650, 191, 21))
        self.progressBarDrawingSpTr.setProperty("value", 0)
        self.progressBarDrawingSpTr.setInvertedAppearance(False)
        self.progressBarDrawingSpTr.setObjectName("progressBarDrawingSpTr")
        self.label_26 = QtWidgets.QLabel(self.centralwidget)
        self.label_26.setGeometry(QtCore.QRect(590, 710, 201, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setUnderline(True)
        self.label_26.setFont(font)
        self.label_26.setObjectName("label_26")
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(590, 750, 111, 20))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.K_toplivo_out = QtWidgets.QLineEdit(self.centralwidget)
        self.K_toplivo_out.setGeometry(QtCore.QRect(700, 750, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.K_toplivo_out.setFont(font)
        self.K_toplivo_out.setReadOnly(True)
        self.K_toplivo_out.setObjectName("K_toplivo_out")
        self.K_stop_engine = QtWidgets.QLineEdit(self.centralwidget)
        self.K_stop_engine.setGeometry(QtCore.QRect(1240, 520, 91, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.K_stop_engine.setFont(font)
        self.K_stop_engine.setReadOnly(False)
        self.K_stop_engine.setObjectName("K_stop_engine")
        self.label_12 = QtWidgets.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(1190, 520, 51, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(320, 650, 41, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.K_step_model = QtWidgets.QLineEdit(self.centralwidget)
        self.K_step_model.setGeometry(QtCore.QRect(370, 650, 51, 22))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.K_step_model.setFont(font)
        self.K_step_model.setReadOnly(False)
        self.K_step_model.setObjectName("K_step_model")
        self.ClearGraph = QtWidgets.QPushButton(self.centralwidget)
        self.ClearGraph.setGeometry(QtCore.QRect(680, 860, 121, 21))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.ClearGraph.setFont(font)
        self.ClearGraph.setStyleSheet("background-color: rgb(233, 122, 127);")
        self.ClearGraph.setObjectName("ClearGraph")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1400, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label.setText(_translate("MainWindow", "Этюды"))
        self.label_23.setText(_translate("MainWindow", "Vz="))
        self.Color.setText(_translate("MainWindow", "color="))
        self.saveDataButton.setText(_translate("MainWindow", "Сохранить изменения"))
        self.z0_field.setText(_translate("MainWindow", "0"))
        self.R_field.setText(_translate("MainWindow", "1"))
        self.label_3.setText(_translate("MainWindow", "R="))
        self.x0_field.setText(_translate("MainWindow", "0"))
        self.label_22.setText(_translate("MainWindow", "Фи:"))
        self.Vy0_field.setText(_translate("MainWindow", "0"))
        self.F_dv_field.setText(_translate("MainWindow", "10"))
        self.y0_field.setText(_translate("MainWindow", "0"))
        self.Phi_field.setText(_translate("MainWindow", "0"))
        self.Color_field.setText(_translate("MainWindow", "1, 1, 0"))
        self.Color_2.setText(_translate("MainWindow", "P="))
        self.label_17.setText(_translate("MainWindow", "y="))
        self.P_field.setText(_translate("MainWindow", "1"))
        self.label_20.setText(_translate("MainWindow", "Vx="))
        self.label_18.setText(_translate("MainWindow", "F_дв:"))
        self.Flamecolor_field.setText(_translate("MainWindow", "1, 0.7, 0.5"))
        self.Vx0_field.setText(_translate("MainWindow", "0"))
        self.label_5.setText(_translate("MainWindow", "m="))
        self.label_21.setText(_translate("MainWindow", "Vy="))
        self.label_16.setText(_translate("MainWindow", "x="))
        self.Color_3.setText(_translate("MainWindow", "Fl_color:"))
        self.Vz0_field.setText(_translate("MainWindow", "0"))
        self.label_19.setText(_translate("MainWindow", "z="))
        self.M_field.setText(_translate("MainWindow", "1"))
        self.label_6.setText(_translate("MainWindow", "Объекты"))
        self.label_7.setText(_translate("MainWindow", "Описание объекта"))
        self.label_8.setText(_translate("MainWindow", "Симуляция"))
        self.label_4.setText(_translate("MainWindow", "Шаг по времени:"))
        self.TStep_field.setText(_translate("MainWindow", "0.01"))
        self.StopButton.setText(_translate("MainWindow", "Стоп"))
        self.K_field.setText(_translate("MainWindow", "10"))
        self.label_11.setText(_translate("MainWindow", "Масштаб:"))
        self.StartButton.setText(_translate("MainWindow", "Поехали!"))
        self.label_24.setText(_translate("MainWindow", "Начальные данные"))
        self.label_25.setText(_translate("MainWindow", "Создание этюда"))
        self.deleteLabel.setText(_translate("MainWindow", "Удаление этюда"))
        self.label_9.setText(_translate("MainWindow", "Название"))
        self.createEtudeButton.setText(_translate("MainWindow", "Создать"))
        self.deleteEtudeButton.setText(_translate("MainWindow", "Удалить выбранный этюд"))
        self.label_28.setText(_translate("MainWindow", "Создание объекта"))
        self.label_29.setText(_translate("MainWindow", "Удаление объекта"))
        self.deleteObjectButton.setText(_translate("MainWindow", "Удалить выбранный объект"))
        self.label_30.setText(_translate("MainWindow", "Параметры"))
        self.label_31.setText(_translate("MainWindow", "Начальные данные"))
        self.createObjectButton.setText(_translate("MainWindow", "Создать"))
        self.label_15.setText(_translate("MainWindow", "Тип:"))
        self.label_32.setText(_translate("MainWindow", "Название:"))
        self.label_2.setText(_translate("MainWindow", "Выбран:"))
        self.label_34.setText(_translate("MainWindow", "Выбран:"))
        self.DrawSpaceshipTrajectory.setText(_translate("MainWindow", "Отрисовка траектории полеты ракеты"))
        self.label_26.setText(_translate("MainWindow", "Вычисленные данные"))
        self.label_10.setText(_translate("MainWindow", "Сожжено топлива:"))
        self.K_toplivo_out.setText(_translate("MainWindow", "0"))
        self.K_stop_engine.setText(_translate("MainWindow", "3500"))
        self.label_12.setText(_translate("MainWindow", "Off F_дв"))
        self.label_13.setText(_translate("MainWindow", "Шагов:"))
        self.K_step_model.setText(_translate("MainWindow", "15000"))
        self.ClearGraph.setText(_translate("MainWindow", "Очистить график"))
from spacewidget import SpaceWidget


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
