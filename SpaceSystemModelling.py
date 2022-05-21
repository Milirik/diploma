from Models import PlanetSystem, Planet, SpaceShip
import numpy as np
from matplotlib.animation import FuncAnimation
import json
import math


class SpaceSystemModelling:
	def __init__(self):
		pass

	def load_info(self, name_etude):
		self.is_load = True
		name_value = f'data_for_analysis'
		name_etude_ = name_etude.replace(" ", "_")
		with open(f"./data_for_analysis/{name_value}_{name_etude_}", "w") as write_file:
			json.dump([self.moveDataCoordinates], write_file)



	def HereAreWeGo(self, is_draw_only_trajectory=False, shag=0.0):
		self.moveDataCoordinates = {
			'x': [],
			'y': [],
			't': [],
			'V': [],
			'cosA': [],
			'r': [],
			'Vr': [],
			'Vfi': [],
			'w': [],
			'R_earth_spitnik': [],
			'vklObert': [],
			'V_sh * Vmoon': [],
			'Угол между финальной скоростью и расстоянием на Землю': []
		}

		self.is_load = False

		def NewPoints(i):
			global phaseObert, Sdobavka, flagStartEngineDobavka, flagStopEngineDobavka, R_earth_spitnik, kToplivaLeft, maxW, signVfi, flagStartEngineMoon, flagStopEngineMoon, max_cnt, t, OnOffEngine, dt, Side, plSystem, ksi, eta, zeta, Vksi, Veta, Vzeta, Dksi, Deta, Dzeta, DVksi, DVeta, DVzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, Dksi_Sh, Deta_Sh, Dzeta_Sh, DVksi_Sh, DVeta_Sh, DVzeta_Sh, F_dv, Alpha, Beta, K_stop_engine
			t += dt

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

			# Увеличиваем шаг интегрирования при приближении к Луне
			if(len(ksi) > 1):
				rast = np.sqrt((ksi_Sh - ksi[1])**2 + (eta_Sh - eta[1])**2 +(zeta_Sh - zeta[1])**2)
				if(rast < 1):
					dt = 0.001
				else:
					dt = 0.01

			# Вывод шага
			if(i % 500 == 0):
				print('[step, t, R_earth_spitnik] ', i, t, R_earth_spitnik)

		
			# Запись данных в файл для анализа
			try:
				if(name_etude == 'Маневр по Оберту.json'):
					r = [ksi_Sh - ksi[1], eta_Sh - eta[1]]
					Vr = (Vksi_Sh * r[0] + Veta_Sh * r[1]) / np.sqrt(r[0]**2 + r[1]**2)
					Vfi = (Vksi_Sh * r[1] - Veta_Sh * r[0]) / np.sqrt(r[0]**2 + r[1]**2)
					w = (Vfi**2)/np.sqrt(r[0]**2 + r[1]**2)
					cosA = (Vksi_Sh * Vksi[1] + Veta_Sh * Veta[1])/(np.sqrt(Vksi_Sh**2 + Veta_Sh **2)*np.sqrt(Vksi[1]**2 + Veta[1] **2))
					
					self.moveDataCoordinates['cosA'].append(cosA)
					self.moveDataCoordinates['Vr'].append(Vr)
					self.moveDataCoordinates['Vfi'].append(Vfi)
					self.moveDataCoordinates['w'].append(w)

					self.moveDataCoordinates['r'].append(np.sqrt((ksi_Sh - ksi[1])**2 + (eta_Sh - eta[1])**2))

					if(flagStartEngineMoon and not flagStopEngineMoon):
						phaseObert+=1
					self.moveDataCoordinates['vklObert'].append(phaseObert)

					V_sh__Vmoon = np.sqrt(Vksi_Sh**2 + Vksi[1]**2) * np.sqrt(Veta_Sh**2 + Veta[1]**2) * cosA
					self.moveDataCoordinates['V_sh * Vmoon'].append(V_sh__Vmoon)

					V_fin__Rearth = math.atan2(ksi_Sh - ksi[0], eta_Sh - eta[0]) - math.atan2(Vksi_Sh, Veta_Sh)
					self.moveDataCoordinates['Угол между финальной скоростью и расстоянием на Землю'].append(V_fin__Rearth)



				
				R_earth_spitnik = np.sqrt((ksi_Sh - ksi[0])**2 + (eta_Sh - eta[0])**2)
				self.moveDataCoordinates['R_earth_spitnik'].append(R_earth_spitnik)



				self.moveDataCoordinates['t'].append(t)
				self.moveDataCoordinates['x'].append(ksi_Sh)
				self.moveDataCoordinates['y'].append(eta_Sh)
				self.moveDataCoordinates['V'].append(np.sqrt(Vksi_Sh**2 + Veta_Sh**2))



				if(R_earth_spitnik > 31.6 and not self.is_load): #
					self.load_info(name_etude)
					print('[load] Success')
					print('[V_fin__Rearth, cosA]', V_fin__Rearth, cosA)

			except BaseException as e:
				print(f'[load] Error: Ошибка в запоминании данных - {e}')



			# Работа двигателя, если использован метод Оберта
			if(name_etude == 'Маневр по Оберту.json'):
				# Включаем двигатель против Луны, когда скорость меняет свой знак
				if(signVfi != (-1 if Vfi < 0 else 1) and t > 70 and not flagStartEngineMoon):
					print('[oh yes..] ', cosA)
					plSystem.get_move_equations(False, True)
					flagStartEngineMoon = True

				# Выключем двигатель против Луны, когда центробежная сила достигает своего максимума
				#if(maxW > w and t > 81.6 and not flagStopEngineMoon):
				if(round(cosA, 3) > shag and t > 81.6 and not flagStopEngineMoon):
					print('[oh no..]')
					plSystem.get_move_equations(False, False)
					flagStopEngineMoon = True
					# print('[V_sh__Vmoon t]', V_sh__Vmoon, t)
				# elif(maxW < w and t > 81.6 and not flagStopEngineMoon):
				# 	maxW = w
			else:
				# Добавочная скорость которую можно использовать по Оберту
				if (t > 89 and not flagStartEngineDobavka):
					plSystem.get_move_equations(False, False, True, Sdobavka)
					kToplivaLeft += Sdobavka
					flagStartEngineDobavka = True
				elif(kToplivaLeft >= Sdobavka and t > 89 and not flagStopEngineDobavka):
					plSystem.get_move_equations(False, False, False, 0)
					flagStopEngineDobavka = True


			# Включение и выключение двигателя по шагу
			for kk in range(len(OnOffEngine)):
				if(i > OnOffEngine[kk]['start'] and not OnOffEngine[kk]['is_started']):
					plSystem.get_move_equations(True, False)
					OnOffEngine[kk]['is_started'] = True
				elif (i > OnOffEngine[kk]['stop'] and not OnOffEngine[kk]['is_stoped']):
					plSystem.get_move_equations(False, False)
					OnOffEngine[kk]['is_stoped'] = True


			# Отображение потраченного топлива
			if(i <= int(K_stop_engine)):
				kToplivaLeft += dt * F_dv
			elif(flagStartEngineMoon and not flagStopEngineMoon): #or (flagStartEngineDobavka and not flagStopEngineDobavka))
				kToplivaLeft += dt * w
			self.K_toplivo_out.setText(str(round(kToplivaLeft, 5)))


			# Изменение расположения объектов и ререндер
			if(is_draw_only_trajectory):
			    plSystem.replace_system_without_draw(ksi, eta, zeta, Vksi, Veta,Vzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh)
			else:
				plSystem.replace_system(ksi, eta, zeta, Vksi, Veta, Vzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, i, OnOffEngine)
				drPlanets = [planet.DrawedPlanet for planet in plSystem.planets]
				drTraces = [planet.DrawedTrace for planet in plSystem.planets]
				return [plSystem.spaceShip.DrawedSpaceShip] + drTraces + drPlanets + [plSystem.spaceShip.DrawedTraceEngineOn] + [plSystem.spaceShip.DrawedTraceEngineOff] + [plSystem.spaceShip.DrawedSpaceShipFlame]
				# + [plSystem.spaceShip.DrawedTraceAfterMoon] + 
				# + [plSystem.spaceShip.DrawedTrace] + \
				       


		global phaseObert, Sdobavka, flagStartEngineDobavka, flagStopEngineDobavka, R_earth_spitnik, name_etude, kToplivaLeft, maxW, signVfi, flagStartEngineMoon, flagStopEngineMoon, max_cnt, t, OnOffEngine, Side, dt, plSystem, ksi, eta, zeta, Vksi, Veta, Vzeta, Dksi, Deta, Dzeta, DVksi, DVeta, DVzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, Dksi_Sh, Deta_Sh, Dzeta_Sh, DVksi_Sh, DVeta_Sh, DVzeta_Sh, F_dv, Alpha, Beta, K_stop_engine
		
		phaseObert = 0
		flagStartEngineDobavka = False
		flagStopEngineDobavka = False
		Sdobavka = 0.74798

		R_earth_spitnik = 0
		kToplivaLeft = 0
		signVfi = 1
		flagStartEngineMoon = False
		flagStopEngineMoon = False
		maxW = 0
		t = 0
		F_dv = 0
		Alpha = 0
		Beta = 0

		# Параметры системы
		dt = float(self.TStep_field.text()) # Шаг интегрирования
		phi = float(self.moonPhi.text()) # Фаза для Луны
		max_cnt = int(self.K_step_model.text()) # Кол-во шагов интегрирования
		name_etude = self.chosenEtudeLabel.text() # Название этюда
		

		razm = 4.216424392e7 # Для обезразмеривания
		koff = 7.29e-5 # Для обезразмеривания


		plSystem = PlanetSystem([])
		for i in self.fileData:
		    if(i['type'] == 'planet'):
		        ksi_, eta_, zeta_ = [i["x"] / razm, i["y"] / razm, i["z"] / razm]

		        V_ksi, V_eta, V_zeta = [i["Vx"]  / (koff * razm), i["Vy"]  / (koff * razm), i["Vz"]  / (koff * razm)]
		        R =  i["R"] / razm
		        M = i["m"]
		        color = i["color"]

		        if(i["name"] == 'Earth'):
		        	ki = 0.9999999998 # Для обезразмеривания
		        else:
		        	# Key
		        	ki = 0.01232376679 # Для обезразмеривания

		        	if(name_etude not in ('Уход с орбиты.json')):
			        	# Меняем фазу Луны
			        	ksi_1 = ksi_ * np.cos(phi) - eta_ * np.sin(phi)
			        	eta_1 = ksi_ * np.sin(phi) + eta_ * np.cos(phi)

			        	V_ksi1 = V_ksi * np.cos(phi) - V_eta * np.sin(phi)
			        	V_eta1 = V_ksi * np.sin(phi) + V_eta * np.cos(phi)

			        	ksi_, eta_, V_ksi, V_eta = ksi_1, eta_1, V_ksi1, V_eta1


		        


		        plSystem.add_new_planet(Planet(ksi_, eta_, zeta_, V_ksi, V_eta, V_zeta, ki, M, R, color))
		    else:
		        ksi_, eta_, zeta_ = [i["x"] / razm, i["y"] / razm, i["z"] / razm]
		        V_ksi, V_eta, V_zeta = [i["Vx"]  / (koff * razm), i["Vy"]  / (koff * razm), i["Vz"]  / (koff * razm)]
		        R =  6 * razm / razm
		        M = i["m"]
		        F_dv =  i["F_dv"]
		        K_stop_engine_ = i["K_stop_engine"]

		        plSystem.add_spaceship(SpaceShip(ksi_, eta_, zeta_, V_ksi, V_eta, V_zeta, M, R, color, F_dv, K_stop_engine_))

		OnOffEngine = [
			{'start': 0, 'stop': int(K_stop_engine_), 'is_started': False, 'is_stoped': False},
		]


		# ===================== Просчитываем траектории полета и получаем вектора ======================= #
		if((len(plSystem.planets) > 0 and hasattr(plSystem, "spaceShip")) or True): # Убрать TRUE
		    plSystem.get_move_equations(True, False)
		    ksi, eta, zeta, Vksi, Veta,Vzeta = plSystem.get_state_vectors()
		    ksi_Sh = plSystem.spaceShip.ksi
		    eta_Sh = plSystem.spaceShip.eta
		    zeta_Sh = plSystem.spaceShip.zeta
		    Vksi_Sh = plSystem.spaceShip.Vksi
		    Veta_Sh = plSystem.spaceShip.Veta
		    Vzeta_Sh = plSystem.spaceShip.Vzeta
		    K_stop_engine = plSystem.spaceShip.K_stop_engine


		# ====================================== Отрисовка графика ====================================== #

		if(is_draw_only_trajectory):
			self.progressBarDrawingSpTr.setValue(0) 
			self.DrawSpaceshipTrajectory.setEnabled(False)

			dt = 0.01
			cnt = 0
			while R_earth_spitnik <= 31.6: # cnt != max_cnt
			    NewPoints(cnt)
			    cnt+=1
			    self.progressBarDrawingSpTr.setValue(cnt*100/max_cnt) 

			self.DrawSpaceshipTrajectory.setEnabled(True)


		Side = float(self.K_field.text()) # Сторона графика. С помощью нее можно увеличить графики
		# self.SpWidget.canvas.axes.set(xlim=[-4.6, -3.7], ylim=[4.8, 3.9], zlim=[-Side, Side])
		self.SpWidget.canvas.axes.set(xlim=[-Side, Side], ylim=[-Side, Side], zlim=[-Side, Side])
		self.SpWidget.canvas.axes.set_title('Это космос')
		self.SpWidget.canvas.axes.set_xlabel('X')
		self.SpWidget.canvas.axes.set_ylabel('Y')
		self.SpWidget.canvas.axes.set_zlabel('Z')

		if(is_draw_only_trajectory):
			TraceKSI_ = []
			TraceETA_ = []
			TraceZETA_ = []

			TraceKSI_2 = []
			TraceETA_2 = []
			TraceZETA_2 = []

			TraceKSI_m = []
			TraceETA_m = []
			TraceZETA_m = []

			TraceKSI_2m = []
			TraceETA_2m = []
			TraceZETA_2m = []

			TraceKSI_.extend(plSystem.spaceShip.TraceKSI[OnOffEngine[0]['start']:OnOffEngine[0]['stop']])
			TraceETA_.extend(plSystem.spaceShip.TraceETA[OnOffEngine[0]['start']:OnOffEngine[0]['stop']])
			TraceZETA_.extend(plSystem.spaceShip.TraceZETA[OnOffEngine[0]['start']:OnOffEngine[0]['stop']])

			TraceKSI_2.extend(plSystem.spaceShip.TraceKSI[OnOffEngine[0]['stop']:])
			TraceETA_2.extend(plSystem.spaceShip.TraceETA[OnOffEngine[0]['stop']:])
			TraceZETA_2.extend(plSystem.spaceShip.TraceZETA[OnOffEngine[0]['stop']:])

			# TraceKSI_m.extend(plSystem.spaceShip.TraceKSI[OnOffEngine[1]['start']:OnOffEngine[1]['stop']])
			# TraceETA_m.extend(plSystem.spaceShip.TraceETA[OnOffEngine[1]['start']:OnOffEngine[1]['stop']])
			# TraceZETA_m.extend(plSystem.spaceShip.TraceZETA[OnOffEngine[1]['start']:OnOffEngine[1]['stop']])

			# TraceKSI_2.extend(plSystem.spaceShip.TraceKSI[OnOffEngine[0]['stop']:OnOffEngine[1]['start']])
			# TraceETA_2.extend(plSystem.spaceShip.TraceETA[OnOffEngine[0]['stop']:OnOffEngine[1]['start']])
			# TraceZETA_2.extend(plSystem.spaceShip.TraceZETA[OnOffEngine[0]['stop']:OnOffEngine[1]['start']])

			# TraceKSI_2m.extend(plSystem.spaceShip.TraceKSI[OnOffEngine[1]['stop']:])
			# TraceETA_2m.extend(plSystem.spaceShip.TraceETA[OnOffEngine[1]['stop']:])
			# TraceZETA_2m.extend(plSystem.spaceShip.TraceZETA[OnOffEngine[1]['stop']:])

			self.SpWidget.canvas.axes.plot(TraceKSI_, TraceETA_, TraceZETA_, ':', color='red')
			self.SpWidget.canvas.axes.plot(TraceKSI_2, TraceETA_2, TraceZETA_2, ':', color='blue')

			# self.SpWidget.canvas.axes.plot(TraceKSI_m, TraceETA_m, TraceZETA_m, ':', color='red')
			# self.SpWidget.canvas.axes.plot(TraceKSI_2m, TraceETA_2m, TraceZETA_2m, ':', color='blue')

		else:
			plSystem.draw(self.SpWidget.canvas.axes)


		self.SpWidget.canvas.show()
		fig = self.SpWidget.canvas.figure

		if(is_draw_only_trajectory):
			pass
		else:
			self.animation = FuncAnimation(fig, NewPoints, interval=dt * 1000, blit=True)


		self.SpWidget.canvas.draw()
