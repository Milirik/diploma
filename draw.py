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



print('[V final empty]', fileDataEmpty[0]['V'][-1])
print('[V final without engine]', fileDataWithoutObert[0]['V'][-1])
print('[V final with engine]', fileDataWithObert[0]['V'][-1])


ax.plot(fileDataEmpty[0]['t'], fileDataEmpty[0]['V'], label='V_empty')
ax.plot(fileDataWithoutObert[0]['t'], fileDataWithoutObert[0]['V'], label='V_without_Obert')
ax.plot(fileDataWithObert[0]['t'], fileDataWithObert[0]['V'], label='V_with_Obert')


# ax.plot(fileData[0]['x'], fileData[0]['y'], marker='o')
# ax.plot(fileData2[0]['x'], fileData2[0]['y'], marker='o', color='red')


# ax.plot(fileData[0]['t'], fileData[0]['Vr'], label='Vr')
# ax.plot(fileData[0]['t'], fileData[0]['Vfi'], label='Vfi')
# ax.plot(fileData[0]['t'], fileData[0]['w'], label='w')

# ax.plot(fileData2[0]['t'], fileData2[0]['Vr'], label='Vr')
# ax.plot(fileData2[0]['t'], fileData2[0]['Vfi'], label='Vfi')
# ax.plot(fileData2[0]['t'], fileData2[0]['w'], label='w')





# ax.plot(fileData[0]['t'], fileData[0]['Vmoon'])
# ax.plot(fileData[0]['t'], fileData[0]['r'])


ax.legend()
plt.show()




# [r min] 0.08378005118286679 -0.4
# [V final with engine] 0.4660889441279318
# [V final without engine] 0.4567847087812945


# [r min] 0.08287299240955395 -0.1
# [V final with engine] 0.49203890886765556
# [V final without engine] 0.4567847087812945


# [r min] 0.08287299240955395 0.0
# [V final with engine] 0.4935340025409224
# [V final without engine] 0.4567847087812945


# [r min] 0.08287299240955395 0.1
# [V final with engine] 0.4899043025149786
# [V final without engine] 0.4567847087812945


# [r min] 0.08287299240955395  0.2
# [V final with engine] 0.4796459730163225
# [V final without engine] 0.4567847087812945

# [r min] 0.08287299240955395
# [V final with engine] 0.3848448730141091
# [V final without engine] 0.4567847087812945
# [Finished in 8.0s]

# [r min] 0.08287299240955395
# [V final with engine] 0.49246480724950614
# [V final without engine] 0.4567847087812945
# [Finished in 11.0s]







# 0,08239 - расстояние от лунs до спутника должно быть


# [V] 0.5220608509310894 : 1.054
# [r] 0.0845624455340885



# [V] 0.6413923967505706 : 1.04 
# [r] 0.033361217785672956


# 0.5858107604618253 : 1.045
# 0.0509448160844163

# 0.50727923676934 : 1.05
# 0.06942298380023768


