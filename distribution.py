import numpy as np
import matplotlib.pyplot as plt
import abapy
from scipy import interpolate
from scipy import ndimage

def make_sigma(path, eps, sigma_factor = 1., eps_factor =1.):
  data = abapy.misc.read_file(path)
  eps_exp   = data[0] * eps_factor         # Experimental strain 
  sigma_exp = data[1] * sigma_factor       # Experimental stress in Pa
  f = interpolate.interp1d(eps_exp, sigma_exp)
  sigma = f(eps)
  return sigma

def derivate(x, y):
  dx = x[1] - x[0]
  y_d = np.gradient(y) 
  return y_d / dx
  

fig = plt.figure(0)
plt.clf()  

E_cu = 210.e9
k_cu = 0.0
path_cu = "S235_sigma_vs_eps.txt"
eps_cu = np.linspace(0.0003, 0.3, 2000)
sigma_cu = make_sigma(path_cu, eps_cu, sigma_factor = 1.e6, eps_factor = 1.)
sigma_d_cu = derivate(eps_cu, sigma_cu)
sigma_d2_cu = derivate(eps_cu, sigma_d_cu)
F_cu = ndimage.gaussian_filter( (E_cu -sigma_d_cu)/(E_cu-k_cu), .1)
p_cu =  derivate(eps_cu, F_cu)

if True:
  ax = fig.add_subplot(3,1,1)
  plt.plot(eps_cu, sigma_cu, label = "Cu")
  plt.grid()
  ax = fig.add_subplot(3,1,2)
  plt.plot(eps_cu, F_cu)
  plt.grid()
  ax = fig.add_subplot(3,1,3)
  plt.plot(eps_cu, p_cu)


E = 70.e9
k = 0.0
path = "Courbe_ref_alu_transv.txt"
eps = np.linspace(0.0005, 0.035, 100)
sigma = make_sigma(path, eps, sigma_factor = 1.e6)
sigma_d = derivate(eps, sigma)
sigma_d2 = derivate(eps, sigma_d)

if True:
  ax = fig.add_subplot(3,1,1)
  plt.plot(eps, sigma, label = "Al")
  plt.legend()
  ax = fig.add_subplot(3,1,2)
  plt.ylabel("$F(\epsilon)$")
  plt.plot(eps, (E -sigma_d)/(E-k))
  ax = fig.add_subplot(3,1,3)
  plt.plot(eps, sigma_d2 / (k-E))



plt.ylabel("$p(\epsilon)$")
plt.grid()
plt.show()


