#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 00:32:19 2023

@author: lyshen
"""

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pymatgen.core import Structure
import pymatgen.analysis.diffusion.analyzer as df
from pymatgen.io.vasp import Xdatcar, Vasprun

time_start = time.time()

#%% Name of file containing data
temp=1200
file_name="vasprun-Br1.xml" 
#%% Load the trajectories of the atom of interest [all functions in this section are pymatgen function]
species_of_interest=('Li') 
number_of_frames_to_discard=10000#this is the number of frames to discard (considered to be of the equilibriation period).  #前10ps不做分析，但我觉得也应该分析才对
vasprun=Vasprun(file_name) #load the vasprun file
structures=vasprun.structures #extract the structures of the MD frames #用Vasprun函数，从vasprun.xml中提取结构
analyzer= df.DiffusionAnalyzer.from_structures(structures[number_of_frames_to_discard:], smoothed=False, specie=(species_of_interest),temperature=temp,time_step=5, step_skip=5,avg_nsteps=100)  #time_step=1, step_skip=1, avg_nsteps=100
#analyzer从10000fs开始，对特定种类原子在特定温度下每隔5fs取样一次；smoothed控制是否平滑MSD，以什么模式平滑处理；
#avg_nsteps确定要平均的时间步数以获得每个时间步的msd，默认值1000通常很好，但一般要求smoothed=”constant”
trajectories = [s.frac_coords for s in analyzer.get_drift_corrected_structures()] #i for i in xxx这样的循环是快速生成具有特定规律的列表
#使用get_drift_corrected_structures函数是为了减少内存使用，因为MD结果有很多结构信息，不需要同时使用
trajectories = np.array(trajectories) #put the trajectories in array format

#%% Load the crystal and MD structures information
#获取Li3OCl的结构，因为所有的空位都被占据，可以确定下Li会占据的位置
crys_site_structure=Structure.from_file('CONTCAR-Livac')   #load the structure that contains all Li sites "defect-free" [Structure.from_file is a pymatgen function] 
MD_initial_structure=vasprun.structures[0] #load an initial MD structure, just to extract indices as below [vasprun.structures are pymatgen functions]

cry_indices= [j for j, site in enumerate(crys_site_structure) if #get indices of Li atoms from the defect-free structure从Li3OCl结构中获取Li的indices
                   site.specie.symbol in species_of_interest]
MD_indices = [j for j, site in enumerate(MD_initial_structure) if #从掺杂结构中获取Li的indices
                   site.specie.symbol in species_of_interest]

#%%获取Li3OCl的分数坐标
lattice = crys_site_structure.lattice #obtain the lattice parameters [pymatgen function]
coor_ref = crys_site_structure.frac_coords[cry_indices]  #get the fractional coordinates of Li sites in the defect free structure 获取Li3OCl中的Li分数坐标
coor_ref = np.array(coor_ref)  #put the coordinates in an array format

#%%估计每个MD帧上的Li原子与Li3OCl中的固定Li位点之间的距离并判断是否是里固定位点最近的Li
dist_mat=np.zeros((len(MD_indices),len(trajectories)))
check_list=[]
for it in range(len(trajectories)):
    dist_matrix = lattice.get_all_distances(coor_ref,  #估计每个 MD 帧上的 Li 原子与无缺陷结构中的固定 Li 位点之间的距离
                                            trajectories[it, MD_indices, :])  #use the lattice.get_all_distances [pymatgen function] to estimate distances between Li atoms at each MD frame and the fixed Li sites in the defect-free structure
    check = dist_matrix == np.min(dist_matrix, axis=0)[None, :] #获取与某位点最近邻Li原子的序号
    #this check will generate a logical matrix. It tells for each site what atom label(s) has a minimal distance with that site
    check_list.append(check) #store this logical matrix for each time instant, in a list.
check_list_mat=np.asarray(check_list) #该矩阵为30000x81x54的尺寸，54为掺杂体系中的Li编号  transform the list to a numpy array for easiness later.
#该矩阵和check_list是一样的，已经存储好了30ps内每1fs的时候距离81个Li3OCl的Li位点最近的Li原子编号

#%%统计check_list_mat中81个Li3OCl的Li位点在30ps内的true数量
n = np.zeros(len(cry_indices))
for i in range(len(check_list_mat)):
    if i != 0:
        for j in range(len(cry_indices)):
            for k in range(len(MD_indices)):
                if check_list_mat[i, j, k] == True and check_list_mat[i, j, k] != check_list_mat[i-1, j, k]:
                    n[j] = n[j] + 1
        
#将Li3OCl中的Li分数坐标和跳跃次数写进txt文件，便于画图
with open('Li_mobility.txt', 'w') as f: # w 为覆盖写入，a 为追加写入
    f.write('x ' + 'y ' + 'z ' + 'num' + '\n')
    for i in range(len(crys_site_structure)):
        f.write(str(coor_ref[i][0]) +' ' + str(coor_ref[i][1]) + ' ' + str(coor_ref[i][2]) +' ' + str(n[i]) + '\n')
f.close()

df = pd.read_table('Li_mobility.txt')
df = pd.DataFrame((x.split(' ') for x in df['x y z num']), columns=['x', 'y', 'z', 'num']) #遇到空格分列
df.to_excel('Li_mobility.xlsx', 'Sheet1', index=False)

#%%绘制三维气泡图
data = pd.read_excel('Li_mobility.xlsx')
x = pd.DataFrame(data, columns=['x'])
y = pd.DataFrame(data, columns=['y'])
z = pd.DataFrame(data, columns=['z'])
num = pd.DataFrame(data, columns=['num'])

x2 = 0.34196937
y2 = 0.34196937
z2 = 0.33785860

fig = plt.figure(figsize=(12, 10))

ax = Axes3D(fig)
ax.scatter(x, y, z, s = num/40, marker = 'o', c = 'red')
ax.scatter(x2, y2, z2, s = 1500, marker = '*', c = 'brown')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#设置白色背景
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

ax.grid(None)
plt.show()

#%%
time_end = time.time()
print('time cost: ', time_end-time_start, ' s')


