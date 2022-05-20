import os
import json
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=[15, 9])
ax = fig.add_subplot(1, 1, 1)


with open(os.path.abspath('.\\data_for_analysis\\' + 'data_for_analysis_Уход_с_орбиты.json') , "r") as read_file1:
    fileDataEmpty = json.load(read_file1)

with open(os.path.abspath('.\\data_for_analysis\\' + 'data_for_analysis_Гравитационный_маневр.json') , "r") as read_file2:
    fileDataWithoutObert = json.load(read_file2)

with open(os.path.abspath('.\\data_for_analysis\\' + 'data_for_analysis_Маневр_по_Оберту.json') , "r") as read_file3:
    fileDataWithObert = json.load(read_file3)



# minT = min(min(fileDataEmpty[0]['r']), min(fileDataWithoutObert[0]['r']), min(fileDataWithObert[0]['r']))

# print(fileData)

#  0.021453565019208812


# print('[r min]', min(fileDataEmpty[0]['r']))
# print('[r min]', min(fileDataWithoutObert[0]['r']))
# print('[r min]', min(fileDataWithObert[0]['r']))







# print('[V final empty]', fileDataEmpty[0]['V'][-1])
# print('[V final without engine]', fileDataWithoutObert[0]['V'][-1])
print('[V final with engine]', fileDataWithObert[0]['V'][-1])

# ax.plot(fileDataEmpty[0]['t'], fileDataEmpty[0]['V'], label='V_empty')
# ax.plot(fileDataWithoutObert[0]['t'], fileDataWithoutObert[0]['V'], label='V_without_Obert')
# ax.plot(fileDataWithObert[0]['t'], fileDataWithObert[0]['V'], label='V_with_Obert')



# ax.plot(fileDataWithObert[0]['t'], fileDataWithObert[0]['V'], label='V_with_Obert')

# ax.plot(fileDataWithObert[0]['vklObert'], fileDataWithObert[0]['V'], label='V_with_Obert')







# ax.plot(fileData[0]['x'], fileData[0]['y'], marker='o')
# ax.plot(fileData2[0]['x'], fileData2[0]['y'], marker='o', color='red')


# ax.plot(fileData[0]['t'], fileData[0]['Vr'], label='Vr')
# ax.plot(fileData[0]['t'], fileData[0]['Vfi'], label='Vfi')
# ax.plot(fileData[0]['t'], fileData[0]['w'], label='w')

# ax.plot(fileData2[0]['t'], fileData2[0]['Vr'], label='Vr')
# ax.plot(fileData2[0]['t'], fileData2[0]['Vfi'], label='Vfi')
# ax.plot(fileData2[0]['t'], fileData2[0]['w'], label='w')


# ax.plot(fileData[0]['t'], fileData[0]['Vmoon'])


# ax.plot(fileData[0]['t'], fileData[0]['Vmoon'])
# ax.plot(fileData[0]['t'], fileData[0]['r'])
ax.plot(
	[
		-1.0, -0.9, -0.8, -0.7, -0.6, -0.5,  -0.4, -0.3, -0.2, -0.1, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
	], 
	[
		0.4382085677096251, 0.4382085677096251, 0.4382085677096251, 0.4382085677096251, 0.4382085677096251, 0.4382085677096251,
		0.4474007737929167, 0.4591624747961364, 0.4685513130241974, 0.4743534002940733, 0.4758975076830916, 0.4721177151887614, 
		0.4614779341291713, 0.4416657177567371, 0.4094633091673678, 0.36200077121897595, 0.30084307759518136, 0.23175058306311605,
		0.1574854899714972
	]
)

ax.legend()
plt.show()




# shag -  -1.0
# [load] Success
# [Vfin, shag] 0.4382085677096251
# shag -  -0.9
# [load] Success
# [Vfin, shag] 0.4382085677096251
# shag -  -0.8
# [load] Success
# [Vfin, shag] 0.4382085677096251
# shag -  -0.7
# [load] Success
# [Vfin, shag] 0.4382085677096251
# shag -  -0.6
# [load] Success
# [Vfin, shag] 0.4382085677096251
# shag -  -0.5
# [load] Success
# [Vfin, shag] 0.4382085677096251
# shag -  -0.4
# [load] Success
# [Vfin, shag] 0.4474007737929167
# shag -  -0.3
# [load] Success
# [Vfin, shag] 0.4591624747961364
# shag -  -0.2
# [load] Success
# [Vfin, shag] 0.4685513130241974
# shag -  -0.1
# [load] Success
# [Vfin, shag] 0.4743534002940733
# shag -  0.0
# [load] Success
# [Vfin, shag] 0.4758975076830916
# shag -  0.1
# [load] Success
# [Vfin, shag] 0.4721177151887614
# shag -  0.2
# [load] Success
# [Vfin, shag] 0.4614779341291713
# shag -  0.3
# [load] Success
# [Vfin, shag] 0.4416657177567371
# shag -  0.4
# [load] Success
# [Vfin, shag] 0.4094633091673678
# shag -  0.5
# [load] Success
# [Vfin, shag] 0.36200077121897595
# shag -  0.6
# [load] Success
# [Vfin, shag] 0.30084307759518136
# shag -  0.7
# [load] Success
# [Vfin, shag] 0.23175058306311605
# shag -  0.8
# [load] Success
# [Vfin, shag] 0.1574854899714972
# shag -  0.9
