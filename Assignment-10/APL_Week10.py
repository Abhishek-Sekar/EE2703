#************************Experiment 10*************************
#*********************Done by Abhishek Sekar EE18B067********



# importing the necessary libraries

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

#Question 1

#reading the file h.csv and extracting the filter coefficients

h_n = np.genfromtxt('h.csv') #filter coefficients

#Question 2

plt.stem(h_n)
plt.title('FIR Filter h[n]')
plt.xlabel(' <- n ->')
plt.ylabel('h[n]')
plt.grid()
plt.show()

#magnitude and phase response
w,H_w = sig.freqz(h_n,1)

plt.subplot(2,1,1)
plt.plot(w,np.abs(H_w))
plt.grid()
plt.ylabel('|H(w)|')
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False) #to remove x ticks

plt.subplot(2,1,2)
plt.plot(w,np.angle(H_w))
plt.grid()
plt.ylabel('phase(H(w))')
plt.xlabel('<- w ->')
plt.suptitle('Magnitude and Phase response of H(w)')
plt.show()

#Question 3
n = np.linspace(1, 2**10, 2**10)
x_n = np.cos(0.2*np.pi*n) + np.cos(0.85*np.pi*n)

plt.stem(x_n[:30])
plt.title('x[n]')
plt.xlabel(' <- n ->')
plt.ylabel('x[n]')
plt.grid()
plt.show()


#Question 4
y_n = np.convolve(x_n,h_n)
plt.stem(y_n[:30])
plt.title('y[n]')
plt.xlabel(' <- n ->')
plt.ylabel('y[n]')
plt.grid()
plt.show()

#Question 5
y_circ_n = np.fft.ifft(np.fft.fft(x_n)*np.fft.fft(np.concatenate((h_n, np.zeros(len(x_n) - len(h_n))))))

plt.stem(y_circ_n[:30])
plt.title('y[n] using circular convolution')
plt.xlabel(' <- n ->')
plt.ylabel('y[n]')
plt.grid()
plt.show()

#Question 6
x_pad_n = np.concatenate((x_n, np.zeros(len(h_n)-1)))
h_pad_n = np.concatenate((h_n, np.zeros(len(x_n)-1)))
y_lin_n = np.fft.ifft(np.fft.fft(x_pad_n)*np.fft.fft(h_pad_n))

plt.stem(y_lin_n[:30])
plt.title('y[n] using DFT based linear convolution')
plt.xlabel(' <- n ->')
plt.ylabel('y[n]')
plt.grid()
plt.show()

#Question 7
with open("x1.csv",'r') as f: #open the csv file in read mode
    
    lines = f.readlines()
    
    f.close()

x_1_n = []
for i,line in enumerate(lines):
    line = line[:-1]
    if(i):
        line = line[:-1]
        line += 'j'
        
    x_1_n.append(complex(line))
     
    
x_1_n = np.asarray(x_1_n,dtype = np.complex)
plt.stem(np.abs(x_1_n[:30]))
plt.title('The Zadoff-Chu Sequence')
plt.xlabel(' <- n ->')
plt.ylabel('|x_1[n]|')
plt.grid()
plt.show()

#correlation
y_corr = np.fft.ifftshift(np.correlate(np.roll(x_1_n,5), x_1_n, "full"))
plt.stem(np.abs(y_corr[:30]))
plt.title('The correlated output')
plt.xlabel(' <- n ->')
plt.ylabel('|y[n]|')
plt.grid()
plt.show()

