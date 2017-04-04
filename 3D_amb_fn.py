import numpy as np
from math import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy import signal


fc = 1000000#0000 #centre freq
tau = 0.000001 #pulse length
B = tau**(-1) #bandwidth
fs = 2*B #sampling freq

t = np.linspace(-tau, tau, fs/500) #TIME VECTOR
fd = np.linspace(-fs/2, fs/2, fs/500) #DOPPLER FREQUENCIES

#title things for plots
if fc == 10000000000:
    title_thing = '10 GHz'
    deci = 50
    pulse = np.sin(2*np.pi*fc*t) #sinusoid pulse
else:
    title_thing = '1 MHz'
    deci = 10
    pulse = np.ones(len(t)) #square

PULSE = np.fft.fft(pulse) #fft
MF = -pulse #time reversed sinusoid or mf response


#plot of pulse
plt.plot(t, pulse, 'r')
plt.title('Simple Square Pulse')# of $f_c$ =' + str(title_thing))
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.show()

#ambiguiy function
matrix = np.zeros(shape=(len(fd),int(2*len(t)-1))) #designing size of matrix
for i in range(len(fd)):
    mid = np.exp(1j * 2 * pi * fd[i] * t)
    shift = pulse * mid
    a = np.convolve(MF, shift)
    matrix[i] = np.abs(np.real(a))

#cut throughs
zer_dp = matrix[int((len(fd)/2)-1)] #cut of surface
zer_td = np.zeros(len(fd)) #cut of surface
for i in range(len(fd)):
    extract = matrix[i][3999]
    zer_td[i] = extract

t1 = np.linspace(-tau, tau, int((2 * len(t))-1)) #setting up timeline for plots
  
#cuts
plt.subplot(2, 1, 1)
plt.plot(t1, zer_dp)
plt.title('Zero-Doppler Cut-through ($f_c$ =' + str(title_thing) + ')')
plt.xlabel('Time Delay (s)')
plt.ylabel('Amplitude')
 
plt.subplot(2, 1, 2)
plt.plot(fd, zer_td)
plt.title('Zero-Time Delay Cut-through ($f_c$ =' + str(title_thing) + ')')
plt.xlabel('Doppler Frequency (Hz)')
plt.ylabel('Amplitude')
plt.show()

#plotting the surface
t1 = np.linspace(-tau, tau, int((2 * len(t))-1)) #setting up timeline
x, y = np.meshgrid(t1, fd) #meshgrid for surface
z = matrix

fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
surf = ax.plot_surface(x, y, z, cmap='jet', rstride=int(deci),
                       cstride=int(deci), linewidth=0,
                       antialiased=False, shade=False)
plt.title('Ambiguity function of Simple Sinusoidal Pulse of $f_c$ =' + str(title_thing))
ax.set_xlabel('Time Delay')
ax.set_ylabel('Doppler Shift')
ax.set_zlabel('Amplitude')
plt.show()
