import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import numpy as np
import sys
from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d

args = sys.argv

# 半球状に点を生成する
def make_points(R):
  points = []
  for i in range(1000):
    t = random.random()
    t = - np.arcsin(1-2*t)
    u = random.random() * 2 * np.pi - np.pi

    x = R*np.cos(t) * np.cos(u)
    y = R*np.cos(t) * np.sin(u)
    z = -abs(R*np.sin(t))
    points.append([x, y, z])
  return points

# 距離
def dis(a, b):
  return pow((a[0]-b[0]),2)+pow((a[1]-b[1]),2)+pow((a[2]-b[2]),2)

# 点群を刈り取る 
def cut_points(points, a, R):
  final = []
  while 1:
    if(len(points)==1):
      
      final.append(points[0])
      break
    if(len(points)==0):
      break
    now_point = points.pop(0)
    new_points = []
    for point in points:
      if dis(now_point, point) > a**2:
        new_points.append(point)
    points = new_points
    final.append(now_point)
  return final

  
# 一回分のシミュレータ
def simulate(R, a):
  #半径Rの半球の表面上に点を間隔a以上で配置
  points = make_points(R)
  points = cut_points(points, a, R)
  #一番低い位置の座標と二番目に低い位置の座標
  return np.array(points)



# N = 10000
# #一回1秒だとして16分
R = int(args[1])
A = float(args[2])
# for i in range(N):
#   z = np.sort(simulate(R, A)[:,2])
#   print('{0},{1},{2},{3}'.format(R,A,R+z[0],R+z[1]))


###　画像生成用　模式図　### 
ps = simulate(R, A)
l = int(args[3])

xs = ps[:,0].tolist()
ys = ps[:,1].tolist()
zs = (1*ps[:,2]).tolist()
print(len(xs)*l)
for i in range(len(xs)):
  for j in range(l):
    print('{0} {1} {2}'.format((R+j+1)*xs[i]/R,(R+j+1)*ys[i]/R,(R+j+1)*zs[i]/R))