import csv
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate

# 打开 CSV 文件
with open('../u_final.csv', 'r') as u:
  
  reader_u = csv.reader(u)
  
  data_u = list(reader_u)
#获得x坐标, 变量u
  X = [float(row[0]) for row in data_u[2:]]
  Uex = [float(row[1]) for row in data_u[2:]]
  Uh = [float(row[2]) for row in data_u[2:]]

  X = np.array(X, dtype = np.float64)
  Uex = np.array(Uex, dtype = np.float64)
  Uh = np.array(Uh, dtype = np.float64)

  plt.plot(X, Uex, label='u_exact')
  plt.plot(X, Uh, label='u_simulated')
  plt.xlabel('x')
  plt.ylabel('value')
  plt.title('result of 1d Helmholtz equation with Galerkin(GS)')
  plt.legend()
  plt.savefig('HOM1D_Helmholtz.png',dpi=300)
  plt.show()

  


