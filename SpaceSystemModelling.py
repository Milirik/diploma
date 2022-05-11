from Models import PlanetSystem, Planet, SpaceShip
import numpy as np
from matplotlib.animation import FuncAnimation


class SpaceSystemModelling:
	def HereAreWeGo(self, is_draw_only_trajectory=False):
		def NewPoints(i):
			global traj, flag, t, dt, plSystem, ksi, eta, zeta, Vksi, Veta, Vzeta, Dksi, Deta, Dzeta, DVksi, DVeta, DVzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, Dksi_Sh, Deta_Sh, Dzeta_Sh, DVksi_Sh, DVeta_Sh, DVzeta_Sh, F_dv, Alpha, Beta, K_stop_engine
			
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

			if(i > int(K_stop_engine) and not flag):
			    plSystem.get_move_equations(True)
			    flag = True

			if(i <= int(K_stop_engine)):
			    self.K_toplivo_out.setText(str(int(i*F_dv)))


			if(is_draw_only_trajectory):
			    print(f'[x] ', i, plSystem.spaceShip.ksi, plSystem.spaceShip.eta, plSystem.spaceShip.zeta, Vksi_Sh, Veta_Sh)
			    plSystem.replace_system_without_draw(ksi, eta, zeta, Vksi, Veta,Vzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh)
			    traj.append([plSystem.spaceShip.ksi, plSystem.spaceShip.ksi, plSystem.spaceShip.ksi])
			else:
				# print(f'[x] ', i, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh)
				print(f'[moon] ', i, ksi, eta, Vksi, Veta)
				plSystem.replace_system(ksi, eta, zeta, Vksi, Veta, Vzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh)
				drPlanets = [planet.DrawedPlanet for planet in plSystem.planets]
				drTraces = [planet.DrawedTrace for planet in plSystem.planets]
				return  [plSystem.spaceShip.DrawedSpaceShip]\
				       + drTraces+drPlanets+ [plSystem.spaceShip.DrawedTrace]


		global traj, flag, t, dt, plSystem, ksi, eta, zeta, Vksi, Veta, Vzeta, Dksi, Deta, Dzeta, DVksi, DVeta, DVzeta, ksi_Sh, eta_Sh, zeta_Sh, Vksi_Sh, Veta_Sh, Vzeta_Sh, Dksi_Sh, Deta_Sh, Dzeta_Sh, DVksi_Sh, DVeta_Sh, DVzeta_Sh, F_dv, Alpha, Beta, K_stop_engine
		t = 0.0
		flag = False

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
		        K_stop_engine_ = i["K_stop_engine"]

		        plSystem.add_spaceship(SpaceShip(ksi_, eta_, zeta_, V_ksi, V_eta, V_zeta, M, R, color, F_dv, K_stop_engine_))

		# F_dv = 0 # Если убрать то будет норм двигатель работать для ракеты

		# ===================== Просчитываем траектории полета и получаем вектора ======================= #
		if((len(plSystem.planets) > 0 and hasattr(plSystem, "spaceShip")) or True): # Убрать TRUE
		    plSystem.get_move_equations(False)
		    ksi,eta,zeta, Vksi, Veta,Vzeta = plSystem.get_state_vectors()
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
			max_cnt = int(self.K_step_model.text())
			traj = []
			while cnt != max_cnt:
			    NewPoints(cnt)
			    cnt+=1
			    self.progressBarDrawingSpTr.setValue(cnt*100/max_cnt) 

			self.DrawSpaceshipTrajectory.setEnabled(True)


		Side = float(self.K_field.text()) # Сторона графика. С помощью нее можно увеличить графики
		self.SpWidget.canvas.axes.set(xlim=[-Side, Side], ylim=[-Side, Side], zlim=[-Side, Side])
		self.SpWidget.canvas.axes.set_title('Это космос')
		self.SpWidget.canvas.axes.set_xlabel('X')
		self.SpWidget.canvas.axes.set_ylabel('Y')
		self.SpWidget.canvas.axes.set_zlabel('Z')

		if(is_draw_only_trajectory):
			self.SpWidget.canvas.axes.plot(plSystem.spaceShip.TraceKSI, plSystem.spaceShip.TraceETA, plSystem.spaceShip.TraceZETA, ':')
		else:
			plSystem.draw(self.SpWidget.canvas.axes)


		self.SpWidget.canvas.show()
		fig = self.SpWidget.canvas.figure

		if(is_draw_only_trajectory):
			pass
		else:
			self.animation = FuncAnimation(fig, NewPoints, interval=dt * 1000, blit=True)

		self.SpWidget.canvas.draw()



