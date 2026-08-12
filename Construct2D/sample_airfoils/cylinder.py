import numpy as np
from matplotlib import pyplot as plt
from math import pi, sin, cos

def main():

  # Length and width

  l = 1.0
  w = 1.0

  # Number of points on curve

  numcrv = 101

  # Allocate x and y

  npoint = numcrv
  x = np.zeros((npoint))
  y = np.zeros((npoint))
  
  # Circle points

  rad = w/2.
  aaxis = w/2;
  baxis = w/8;
  cenx = rad
  ceny = 0.
  dtheta = 2*pi/float(numcrv-1)
  theta0 = 0.
  for i in range(0,numcrv):
    #x[i] = cenx + rad*cos(theta0 + i*dtheta)
    #y[i] = ceny + rad*sin(theta0 + i*dtheta)
    x[i] = cenx + aaxis*cos(theta0 + i*dtheta)
    y[i] = ceny + baxis*sin(theta0 + i*dtheta)

  # Plot

  ax = plt.subplot(111)
  ax.set_aspect('equal','datalim')
  ax.plot(x,y,'-x')
  plt.show()

  # Write to file

  f = open('cylinder.dat','w')
  f.write('Cylinder\n')
  for i in range(0,npoint-1):
    f.write(str(x[i]) + '  ' + str(y[i]) + '\n')
  f.close()

if __name__ == "__main__":
  main()
