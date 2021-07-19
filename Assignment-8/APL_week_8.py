#**************************EXPERIMENT 8:DFT****************************************
#**************************By Abhishek Sekar***************************************
from pylab import *  #removes the issue of elaborating library functions but makes code slow

#Function to select different functions
def select(t, n):
    if(n == 1):
        return sin(5*t)      #function in example 1
    elif(n == 2):
        return (1+0.1*cos(t))*cos(10*t)#function in example 2
    elif(n == 3):
        return pow(sin(t), 3)
    elif(n == 4):
        return pow(cos(t), 3)
    elif(n == 5):
        return cos(20*t + 5*cos(t))
    else:
        return exp(-pow(t, 2)/2)   #Gaussian


'''
Function to find Discrete Fourier Transform
 l_lim,u_lim    -> lower & upper limit for time
 points         -> Sampling rate
 n              -> code for function above(1,6)
 norm_factor    ->  only for Gaussian function
                   it is given as parameter, it is the normalizing factor.
'''


def DFT(l_lim, u_lim, points,n,norm_Factor=0):
    t =np.linspace(l_lim, u_lim,points+1)[:-1]
    y =select(t, n)
    N = points

    # DFT for gaussian function
    # norm_factor is multiplying constant to DFT

    if(norm_Factor != 0):
        Y = fftshift((fft(y))*norm_Factor)
    else:
        # normal DFT for periodic functions
        Y = fftshift(fft(y))/(N)

    w_lim = (2*pi*N/((u_lim-l_lim)))  #freq domain
    w = linspace(-(w_lim/2), (w_lim/2), (points+1))[:-1]
    return t, Y, w


'''
Function to plot Magnitude and Phase spectrum for given function
Arguments:

 thresh_val     -> value above which phase is made zero
 xlims,ylims    -> limits for x&y axis for spectrum
 plt_title -> title of plot 
'''


def plot_DFT(t, Y, w, thresh_val,plt_title,xlims=[-15,15],ylims=None):

    subplot(2, 1, 1)
    plt.plot(w, abs(Y), lw=2)
    plt.xlim(xlims)
    if(ylims != None):
        plt.ylim(ylims)
    plt.xlabel("w ->")
    plt.ylabel("|Y(jw)| ->")
    plt.title(plt_title)
    plt.grid()

    ax = subplot(2, 1, 2)
    ii = where(abs(Y) > thresh_val)
    plt.plot(w[ii], angle(Y[ii]), 'go', lw=2)

    if(ylims != None):
        plt.ylim(ylims)

    plt.xlim(xlims)
    plt.ylabel("phase(Y(jw)) ->")
    plt.xlabel("w ->")
    plt.grid()
    plt.show()


'''
DFT for sin(5t) computed in incorrect way
* like without normalizing  factor
* without centering fft of function to zero
'''

x = linspace(0, 2*pi, 128)
y = sin(5*x)
Y = fft(y)
subplot(2, 1, 1)
plt.plot(abs(Y), lw=2)
plt.title("Figure 1 : Incorrect Spectrum of sin(5t)")
plt.ylabel("|Y(jw)| ->")
plt.grid()

subplot(2, 1, 2)
plt.plot(unwrap(angle(Y)), lw=2)
plt.xlabel("w -> ")
plt.ylabel("phase(Y(jw)) ->")
plt.grid()
plt.show()


#proper computation of DFT for Sin(5t)
t, Y, w = DFT(0, 2*pi, 128,1)
plot_DFT(t, Y, w,1e-3,"Figure 2: Spectrum of sin(5t)")

#unstretched t-axis calculation for (1+0.1cos(t))cos(10t)
t, Y, w = DFT(0, 2*pi, 128,2)
plot_DFT(t, Y, w, 1e-4,"Figure 3: Incorrect Spectrum of (1+0.1cos(t))cos(10t)")

#stretched t-axis calculation for refined plot
t, Y, w = DFT(-4*pi, 4*pi, 512,2)
plot_DFT(t, Y, w, 1e-4, "Figure 4 : Spectrum of (1+0.1cos(t))cos(10t)")

#subdivision 2: plot for sin^3(t) 
t, Y, w = DFT(-4*pi, 4*pi, 512,3)
plot_DFT(t, Y, w, 1e-4, "Figure 5: Spectrum of sin^{3}(t)")

#subdivision 2: plot for cos^3(t) 
t, Y, w = DFT(-4*pi, 4*pi, 512,4)
plot_DFT(t, Y, w, 1e-4,"Figure 6: Spectrum of cos^{3}(t)")

#subdivision 3: plot for cos(20t + 5cos(t))
t, Y, w = DFT(-4*pi, 4*pi, 512,5)
Xlims = [-40, 40]
plot_DFT(t, Y, w, 1e-3,"Figure 7: Spectrum of cos(20t+5cos(t))",Xlims)

#Subdivision 4: For the Gaussian
# initial window_size and sampling rate defined
window_size = 2*pi
sampling_rate = 128
# tolerance for error
tol = 1e-15

# normalisation factor derived
norm_factor = (window_size)/(2*pi*(sampling_rate))

#for loop to identify where the minimum error occurs
#window_size and sampling_rate increased to improve accuracy and overcome aliasing problems

for i in range(1, 10):

    t, Y, w = DFT(-window_size/2, window_size/2,
                      sampling_rate,6, norm_factor)

    # actual transform of the gaussian
    actual_Y = (1/sqrt(2*pi))*exp(-pow(w, 2)/2)
    #comparing both to get error
    error = (mean(abs(abs(Y)-actual_Y)))
    
    if(error < tol):
        print("\nAccuracy of the DFT is: %g and Iterations took: %g" %
              ((error, i)))
        print("Best Window_size: %g , Sampling_rate: %g" %
              ((window_size, sampling_rate)))
        break
    else:
        window_size = window_size*2
        sampling_rate = (sampling_rate)*2
        norm_factor = (window_size)/(2*pi*(sampling_rate))


Xlims = [-10, 10]
plot_DFT(t, Y, w, 1e-2,r"Figure 8: Spectrum of $e^{-\frac{t^{2}}{ 2}}$",Xlims)

# Plotting actual DFT of Gaussian
plot(w, abs(actual_Y),
     label=r"$\frac{1}{\sqrt{2}\pi} e^{\frac{\ - \omega ^{2}}{2}}$")
title("Exact Fourier Transform of Gaussian")
xlim([-10, 10])
ylabel("Y(jw) ->")
xlabel("w ->")
grid()
legend()
show()