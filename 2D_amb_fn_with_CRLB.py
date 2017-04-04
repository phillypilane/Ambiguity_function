from math import *
import cmath as cmath
import numpy as np
import matplotlib.pyplot as plt
import random as random
from scipy import signal

B = 1000000 #bandwidth
tau = 0.001 #pulse length
rate = 10000 #sample freq
t = np.linspace(0, tau/2, rate/2) #actual chirp length

#LFM chirp (unpadded) creation
unpadded = []
for i in range (len(t)):
    yyy = cmath.exp(np.complex(0,1) * np.pi * (B/tau) * t[i]**2)
    unpadded.append(yyy)

#plotting the chirp
plt.plot(t, unpadded)
plt.xlabel('Time (s)')
plt.ylabel('Magnitude')
plt.title('LFM Chirp (Not Zero-Padded)')
plt.show()

#Padding of the chirp
pad = np.zeros(rate/4)
half = pad.tolist()
v_tx = half + unpadded + half #adding pad on either side

timeline = np.linspace(0, tau, len(v_tx)) #real timeline
timeline2 = np.linspace(-tau/2, tau/2, len(v_tx)) #+/- timeline for display only

#plotting padded chirp
plt.plot(timeline, v_tx)
plt.xlabel('Time (s)')
plt.ylabel('Magnitude')
plt.title('LFM Chirp (Zero-Padded)')
plt.show() 

#ambiguity function (zero doppler i.e. fd = 0)
zero_dop = []
for i in range(len(t)):
    con = unpadded + half + half #creating chirp that will shift through timeline
    x_star = np.conj(con)
     
    tr = np.roll(x_star, i) #shifting it through time as loop iterates
     
    ty = np.multiply(v_tx,tr)
    zero_dop.append(abs(sum(ty)))
 
#plotting Ambiguity fn of zero doppler
plt.plot(timeline2, (half + zero_dop + half))
plt.xlabel('Time (s)')
plt.ylabel('Magnitude')
plt.title('Ambiguity Function (Zero Doppler)')
plt.show()

# #ambiguity function (zero time delay i.e. t = 0)
doppler = B/2 #doppler frequency limit
fd = np.linspace(-doppler, doppler, rate) #doppler frequency range
 
zero_time = []
for i in range(len(fd)):
    
    num = sin(np.pi * fd[i] * tau)
    den = np.pi * fd[i] * tau
    frac = abs(num/den)
    zero_time.append(frac)

#plotting zero time delay
plt.plot(fd, zero_time)
plt.xlabel('Doppler Frequency (Hz)')
plt.ylabel('Magnitude')
plt.title('Ambiguity Function (Zero Time Delay)')
plt.show()

#Hamming window
ham = np.hamming(len(t)) #calling straight from numpy
hamming = ham.tolist() #array to list
hamming_padded = half + hamming + half #zero padding

#plotting pure hamming window
plt.plot(timeline,hamming_padded)
plt.xlabel('Time (s)')
plt.ylabel('Magnitude')
plt.title('Hamming Window - Zero Padded')
plt.show()

#creating window chirp
H = np.conjugate(v_tx) #complex conjugate of tx signal
windowed_mf = np.multiply(H, hamming_padded) #point multiplication as suggested

#plot the windowed chirp
plt.plot(timeline, windowed_mf)
plt.xlabel('Time (s)')
plt.ylabel('Magnitude')
plt.title('Hamming Windowed Chirp')
plt.show()

#tapering the matched filter response
amb_fn = []
for i in range(len(timeline)):
    con = unpadded + half + half
    x_star = (windowed_mf)
    tp = np.roll(x_star, i -int(len(x_star)/2))
    
    app = np.multiply(v_tx, tp)
    amb_fn.append(abs(sum(app)))

#plot the tapered MF R
plt.plot(timeline, amb_fn)
plt.title('Match Filter Response (Hamming Tapered)')
plt.xlabel('Time (s)')
plt.ylabel('Magnitude')
plt.show()

#Cramer rao lower bound
f0 = np.linspace(0.01,0.49,500) #working between normalised freq
sigma = 1 #sigma or standard dev
phase = 0 #phase
A = 1 #amplitude
N = 15

var = []
for i in range(len(f0)):
    f = f0[i]
    
    den_part = []
    for n in range(N - 1):
        uuu = float(2 * pi * n * sin(2*pi* f * n + phase))**2
        den_part.append(uuu)

    frac = (sigma**2) / ((A**2) * sum(den_part))
    var.append(frac)

#plotting the CRLB
plt.plot(f0, var)
plt.title('Cramer-Rao Lower Bound - Var VS Normalised Frequency')
plt.xlabel('Normalised Frequency')
plt.ylabel('Variance')
plt.show()
