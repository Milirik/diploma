import os
import json
import numpy as np
import matplotlib.pyplot as plt

with open(os.path.abspath('.\\data_for_analysis\\' + 'DataForAnalysis.json') , "r") as read_file:
    fileData = json.load(read_file)

# with open(os.path.abspath('.\\data_for_analysis\\' + 'DataForAnalysis_save.json') , "r") as read_file2:
#     fileData2 = json.load(read_file2)

# print(fileData)


fig = plt.figure(figsize=[15, 9])
ax = fig.add_subplot(1, 1, 1)
# ax.plot(fileData[0]['x'], fileData[0]['y'], marker='o')
# ax.plot(fileData2[0]['x'], fileData2[0]['y'], marker='o', color='red')

ax.plot(fileData[0]['t'], fileData[0]['V'])

ax.plot(fileData[0]['t'], fileData[0]['r'])

print('[V]', fileData[0]['V'][-1])
print('[r]', min(fileData[0]['r']))



plt.show()

# 0,08239


# [V] 0.5220608509310894 : 1.054
# [r] 0.0845624455340885



# [V] 0.6413923967505706 : 1.04 
# [r] 0.033361217785672956


# 0.5858107604618253 : 1.045
# 0.0509448160844163

# 0.50727923676934 : 1.05
# 0.06942298380023768


