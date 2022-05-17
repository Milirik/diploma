import os
import json
import numpy as np

with open(os.path.abspath('.\\data_for_analysis\\' + 'DataForAnalysis.json') , "r") as read_file:
    fileData = json.load(read_file)


max_v = 0
max_i = 0
min_v = 9999999999999999999
min_i = 0


for Vx, Vy in zip(fileData[0]['Vx'], fileData[0]['Vy']):
	V = np.sqrt(Vx['Vx']**2 + Vy['Vy']**2)
	print('[x]', Vx['i'], V)
    
	if max_v < V:
		max_v = V
		max_i = Vx['i']

	if min_v > V:
		min_v = V
		min_i = Vx['i']

print(f'max_v: {max_v}, i: {max_i}')
print(f'min_v: {min_v}, i: {min_i}')

# max_v: 0.514312552259352, i: 5665
# min_v: 0.2225524154423358, i: 13594

	