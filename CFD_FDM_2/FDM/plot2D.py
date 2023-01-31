import csv
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate


def plot_contour(Points,Strain,minX,maxX,minY,maxY,minZ,maxZ,figName):
#生成contour的网格
    X = np.linspace(minX, maxX, 256)
    Y = np.linspace(minY, maxY, 256)
    X1, Y1 = np.meshgrid(X, Y)
#对data进行插值，Points为二维数组,比如2D的情况下的为(：，2)；输入的Data(Strain)为一维数组
    Z = interpolate.griddata(Points, Strain, (X1, Y1), method='linear', fill_value=0)
    Z[np.where (Z < minZ)] = minZ                                                                                                                                                                                                                                    
    Z[np.where (Z > maxZ)] = maxZ

    fig, ax = plt.subplots(figsize=(12, 12))    
    levels = np.linspace(minZ,maxZ,60) 
    cset1 = ax.contourf(X1, Y1, Z, levels,cmap=cm.jet) 
  
    ax.set_title("w_result t=0.1s",size=18)  
    ax.set_xlim(minX, maxX)
    ax.set_ylim(minY, maxY)
    ax.set_xlabel("X(m)",size=18)
    ax.set_ylabel("Y(m)",size=18)

    cbar = fig.colorbar(cset1)
    cbar.set_label('w', size=18)
    cbar.set_ticks(np.linspace(minZ,maxZ,10))

    fig.savefig(figName+".png", bbox_inches='tight',dpi=300,pad_inches=0.1)
    plt.show()  

    return()

def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2


# 打开 CSV 文件
with open('../w_final.csv', 'r') as f:
  
  reader = csv.reader(f)
  
  data = list(reader)
#获得x和y坐标
  X = [float(row[0]) for row in data]
  Y = [float(row[0]) for row in data]
  X = np.array(X, dtype = np.float64)
  Y = np.array(Y, dtype = np.float64)
#形成x与y坐标构成的网格
  X_mesh, Y_mesh = np.meshgrid(X, Y)
  Points = np.column_stack((Y_mesh.flatten(), X_mesh.flatten()))
#获得data，并平铺为一维数组
  data = np.delete(data, [0, -1], axis=1)
  data_array = np.array(data, dtype = np.float64)
  data_array = data_array.flatten()

  dataMax=np.max(data_array)
  dataMin=np.min(data_array)

plot_contour(Points,data_array,min(X),max(X),min(Y),max(Y),dataMin,dataMax,"contour")
