import numpy as np
import matplotlib.pyplot as plt
import abapy
from scipy import interpolate, signal 


def derivate(x, y):
  dx = x[1] - x[0]
  y_d = np.gradient(y) 
  return y_d / dx


class TensileTest(object):
  
  def __init__(self, path, eps_min = None, eps_max = None, Np = None, eps_factor = 1., sigma_factor = 1., E = 1., k = 0.):
    data = abapy.misc.read_file(path)
    eps_raw   = data[0] * eps_factor         # Experimental strain 
    sigma_raw = data[1] * sigma_factor       # Experimental stress in Pa
    if eps_min == None: eps_min = eps_raw.min()
    if eps_max == None: eps_max = eps_raw.max()
    if Np == None: Np = len(eps_raw)
    self.eps = np.linspace(eps_min, eps_max, Np) 
    self.sigma = interpolate.interp1d(eps_raw, sigma_raw)(self.eps)
    self.E = E
    self.k = k
    
  def filter_sigma(self, core = 1.e3):
    eps = self.eps
    deps = eps[1] - eps[0]
    size = int(core  / deps)
    if (size % 2 == 0): size += 1
    self.sigma = signal.medfilt( self.sigma, size)
    
  def get_dsigma(self):    
    return derivate(self.eps, self.sigma)
  dsigma = property(get_dsigma)  
      
  def get_ddsigma(self):    
    return derivate(self.eps, self.dsigma)
  ddsigma = property(get_ddsigma)      

  def get_cdf(self):
    E = self.E
    k = self.k
    dsigma = self.dsigma
    return (E -dsigma)/(E - k)
  cdf = property(get_cdf)  
    
    
  def get_pdf(self):
    E = self.E
    k = self.k
    ddsigma = self.ddsigma
    return ddsigma / (k - E)
  pdf = property(get_pdf)    








fig = plt.figure(0)
plt.clf()  

ax1 = fig.add_subplot(3,1,1)
plt.grid()
plt.xlim(0., 0.1)
plt.ylabel("$\sigma$")
plt.legend()
ax2 = fig.add_subplot(3,1,2)
plt.grid()
plt.legend()
plt.ylabel("CDF")
plt.xlim(0., 0.1)
ax3 = fig.add_subplot(3,1,3)
plt.grid()
plt.legend()
plt.ylabel("PDF")
plt.xlim(0., 0.1)
 
data = {} 
data["Cu"] = TensileTest(path = "Cu_sigma_vs_eps.txt", 
                         eps_min = 1.e-4,
                         eps_max = None,
                         eps_factor = 1.,
                         sigma_factor = 1.,
                         E = 103.e3,
                         k = 0.,
                         Np = 200)

data["Al 6060"] = TensileTest(path = "Courbe_ref_alu_transv.txt", 
                         eps_min = 1.5e-3,
                         eps_max = None,
                         eps_factor = 1.,
                         sigma_factor = 1.,
                         E = 69.5e3,
                         k = 0.,
                         Np = 100)

data["S235 Steel"] = TensileTest(path = "S235_sigma_vs_eps.txt", 
                         eps_min = 1.e-4,
                         eps_max = .3,
                         eps_factor = 1.,
                         sigma_factor = 1.,
                         E = 210e3,
                         k = 0.,
                         Np = 500)

colors = ["r", "b", "g", "m", "y", "c"]
keys = data.keys()
for i in  xrange(len(keys)):
  k  = keys[i]
  tt = data[k]
  #ax1.plot(tt.eps, tt.sigma, color = colors[i], marker = ",", linestyle = "", label = k)
  #tt.filter_sigma(core = .02)
  ax1.plot(tt.eps, tt.sigma, color = colors[i], marker = "", linestyle = "-", label = k)
  ax2.plot(tt.eps, tt.cdf, color = colors[i], marker = "", linestyle = "-", label = k)
  ax3.plot(tt.eps, tt.pdf, color = colors[i], marker = "", linestyle = "-", label = k)
ax1.legend()
ax3.legend()
ax2.legend()
#plt.tight_layout()
plt.savefig("Experimental_distributions.pdf")


