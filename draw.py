import os
import json
import numpy as np
import matplotlib.pyplot as plt

with open(os.path.abspath('.\\data_for_analysis\\' + 'DataForAnalysis.json') , "r") as read_file:
    fileData = json.load(read_file)

with open(os.path.abspath('.\\data_for_analysis\\' + 'DataForAnalysis_save.json') , "r") as read_file2:
    fileData2 = json.load(read_file2)

# print(fileData)

print(fileData[0]['V'][-1])

fig = plt.figure(figsize=[15, 9])
ax = fig.add_subplot(1, 1, 1)
# ax.plot(fileData[0]['x'], fileData[0]['y'], marker='o')
# ax.plot(fileData2[0]['x'], fileData2[0]['y'], marker='o', color='red')

ax.plot(fileData[0]['t'], fileData[0]['V'])


plt.show()




