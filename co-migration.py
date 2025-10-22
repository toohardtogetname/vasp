# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 20:14:01 2025

@author: lshenaf
"""


import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Xdatcar

#dirname = r'D:\OneDrive - HKUST Connect\Attachments\halide\oxyhalide\NbCl5-based\AIMD\0.50Li2O\MLFF_IVDW\\'
with open('co-migration proportion.csv', 'a') as f1:
    f1.write('Li_index' + '|' + 'neighbor atom' + '|' + 'neighbor atom_index' + '|' + 'proportion' + '|' + '\n')#)

# 读取 XDATCAR 文件
xdatcar = Xdatcar('XDATCAR')
traj = xdatcar.structures


# 提取初始结构
initial_struct = traj[0]

# 找到初始结构中的 Li 原子索引（假设元素符号是 Li）
Li_indices = []
for i, atom in enumerate(initial_struct):
    if atom.species_string == 'Li':
        Li_indices.append(i)

'''
# 指定感兴趣的 Li 原子的索引（例如第一个 Li 原子）
Li_index = Li_indices[20]
cutoff = 4  # 截断半径，单位为Å
neighborlist = []
'''

for x in range(len(Li_indices)):
    Li_index = Li_indices[x]
    cutoff = 4  # 截断半径，单位为Å
    neighborlist = []

    for i, atom in enumerate(initial_struct):
        if atom.species_string != 'Li':
            distance = initial_struct.get_distance(Li_index, i)
            if distance <= cutoff:
                print(f"atom_index={i+1}, {atom.species_string}, distance={distance}") #i+1是因为程序从0开始，而vesta从1开始
                neighborlist.append(i)
    
    print(f"The neighbor list of Li{Li_index} with cutoff {cutoff}Å is: ", neighborlist)
    
    for j in range(len(neighborlist)):  
        num = 0
        for k in range(len(traj)):
            distance = traj[k].get_distance(Li_index, neighborlist[j])
            if distance <= cutoff:
                num = num + 1
        proportion = num/len(traj)*100
        print(f"The proportion of Li-{Li_index+1} with atom{neighborlist[j]+1}-{traj[k][neighborlist[j]].species_string} is: ", proportion)
        with open('co-migration proportion.csv', 'a') as f1:
            f1.write(str(Li_index+1) + '|' + traj[k][neighborlist[j]].species_string + '|' + str(neighborlist[j]+1) + '|' + str(proportion) + '%' + '|' + '\n')
        

#最后把所有的Li周围的原子距离及距离保持在4Å以内的时间步的比例统计出来
    

