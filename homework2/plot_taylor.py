import matplotlib.pyplot as plt
import numpy as np

def f0(x):
	return x
def f1(x):
	return x-x**3/6
def f2(x):
	return x-x**3/6+x**5/120
def f3(x):
	return x-x**3/6+x**5/120-x**7/5040
def f4(x):
	return np.sin(x)

x=np.arange(-10.0,10.0,0.01)
plt.plot(x,f0(x),x,f1(x),x,f2(x),x,f3(x),x,f4(x))
plt.axis([-10,10,-10,10])
plt.show()



