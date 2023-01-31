import csv
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate

# 打开 CSV 文件
with open('../rho_final.csv', 'r') as rho, open('../u_final.csv', 'r') as u, open('../e_final.csv', 'r') as e:
  
  reader_rho = csv.reader(rho)
  reader_u = csv.reader(u)
  reader_e = csv.reader(e)
  
  data_rho = list(reader_rho)
  data_u = list(reader_u)
  data_e = list(reader_e)
#获得x坐标, 变量rho, u, e
  X = [float(row[0]) for row in data_rho[1:]]
  Rho = [float(row[11]) for row in data_rho[1:]]
  U = [float(row[11]) for row in data_u[1:]]
  E = [float(row[11]) for row in data_e[1:]]

  X = np.array(X, dtype = np.float64)
  Rho = np.array(Rho, dtype = np.float64)
  U = np.array(U, dtype = np.float64)
  E = np.array(E, dtype = np.float64)

  plt.plot(X, Rho, label='rho')
  plt.plot(X, U, label='u')
  plt.plot(X, E, label='e')
  plt.xlabel('x')
  plt.ylabel('value')
  plt.title('result of Euler equation with HLLC_Riemann_Solver')
  plt.legend()
  plt.savefig('HLLC_Riemann_Solver.png',dpi=300)
  plt.show()

  


