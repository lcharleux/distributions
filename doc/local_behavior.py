import numpy as np
import matplotlib.pyplot as plt
import os

#Latex preamble
from matplotlib import rc
import matplotlib
fig_width = 140 * 0.0393700787 # 90, 140, 190
fig_size = (fig_width, fig_width/1.6)

rc('text', usetex=True)
rc('text.latex', preamble = 
    '''\usepackage{amsmath}
    \usepackage{yfonts}
    \usepackage[T1]{fontenc}
    \usepackage[latin1]{inputenc}
    \usepackage{txfonts}
    \usepackage{lmodern}
    \usepackage{blindtext} 
    \usepackage{wasysym}
    \usepackage[bitstream-charter]{mathdesign}''')
rc('font', family='serif', weight='normal', style='normal')
params = {'backend': 'pdf',
          'axes.labelsize': 12,
          'text.fontsize': 12,
          'legend.fontsize': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'figure.figsize': fig_size,
          'axes.unicode_minus': True}
matplotlib.rcParams.update(params)

E = 1.
n = .25
sy = .25


k = (E**-1 + n**-1)**-1
ey = sy /E

xe = np.linspace(0., ey, 3)
xp = np.linspace(ey, 1., 3)

se = E * xe
sp = sy + k * (xp-sy/E)

fig = plt.figure(0)
plt.clf()
ax = fig.add_subplot(1,1,1)
plt.plot(xe, se, "bo-", label = "Elastic")
plt.plot(xp, sp, "rs-", label = "Elastic-plastic")
plt.grid()
plt.legend(loc = "best")
plt.xticks([0, sy/E],["$0$", "$\epsilon_y = \sigma_y/E$"])
plt.yticks([0, sy],["$0$", "$\sigma_y$"])
plt.xlabel("Total Strain, $\epsilon = \epsilon_e + \epsilon_p$")
plt.ylabel("Stress, $\sigma$")

ax.annotate('$\sigma = E \epsilon$', xy=(ey/2., sy/2.), xytext=(ey/2.+.3, sy/2.-.1), ha='center', va = 'center',
            arrowprops=dict(arrowstyle="<|-",
                                shrinkA=5,
                                shrinkB=5,
                                fc="w", ec="b",
                                connectionstyle="angle3,angleA=0,angleB=-45",
                                ),)

xxp = (ey +1.)/2.
yyp = sy + (xxp - ey)*k                                
ax.annotate('$\sigma = \sigma_y + k (\epsilon -\epsilon_y)$', xy=(xxp, yyp), xytext=(xxp +.1, yyp-.1), ha='center', va = 'center',
            arrowprops=dict(arrowstyle="<|-",
                                shrinkA=5,
                                shrinkB=5,
                                fc="w", ec="r",
                                connectionstyle="angle3,angleA=0,angleB=-45",
                                ),)                                
plt.tight_layout()
plt.savefig("local_behavior.pdf")  
