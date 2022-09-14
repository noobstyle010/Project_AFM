from sys import stdin
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import sys 

args = sys.argv
File_name = args[1]

def zscore(x, axis = None):
  xmean = x.mean(axis=axis, keepdims=True)
  xstd  = np.std(x, axis=axis, keepdims=True)
  zscore = (x-xmean)/xstd
  return zscore

input = stdin.readline
###  具体的な処理　###
sn,pn = map(int, input().split())
nx,ny,nz = map(int, input().split())
dx,dy,dz = map(float, input().split())
F = np.ones((nx,ny))
ims = []
for i in range(sn):
  x_,y_,z_ = map(float, input().split())
for i in range(pn):
  x_,y_,z_ = map(float, input().split())
num = 0
for iy in range(ny):
    for ix in range(nx):
        x_,y_,z_,f_ =  map(float, input().split())
        F[ix][iy] = f_
fig = plt.figure(figsize=(0.5,0.5))
f = F
plt.axis("off")
norm =  mcolors.TwoSlopeNorm(vmin=-1*abs(max(f.min(),f.max(),key=abs)),vcenter=0., vmax=abs(max(f.min(),f.max(),key=abs)))
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
plt.imshow(F,norm=norm,cmap='gray')
plt.savefig(args[1] + ".jpg")


