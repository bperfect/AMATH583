import matplotlib
matplotlib.use('Agg')
from matplotlib import image, animation
from matplotlib import pyplot as plt
import JSAnimation_frametools as J
from JSAnimation import IPython_display
import numpy as np
import glob


plotdir='_plots'
J.make_plotdir(plotdir, clobber='true')
fid = open('input_data.txt')
n=int(fid.readline())
tfinal=float(fid.readline())
fid.readline()
nsteps=int(fid.readline())

ax=plt.axes(xlim=(0,np.pi),ylim=(-1,1))
line, = ax.plot([],[],lw=2)
line.set_data(range(n),range(n))

data = open('solution.txt')
for i in range(nsteps+1):
    soln = []
    plt.clf()
    for j in range(n+2):
        soln.append(float(data.readline()))
    plt.axes(xlim=(0,np.pi),ylim=(-1,1))
    plt.plot(np.arange(0,np.pi,np.pi/(n+2)),soln)
    plt.plot(np.arange(0,np.pi,np.pi/(n+2)),soln,'bo')
    t="%.3f" %(i*tfinal/nsteps)
    plt.title('time = '+t)
    J.save_frame(i)
J.make_html(J.make_anim(plotdir))
    
    

