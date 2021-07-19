#************************Experiment 9*************************
#*********************Done by Abhishek Sekar EE18B067*********


from pylab import *
from mpl_toolkits.mplot3d import Axes3D

#worked examples
#sin(sqrt(2)t) attempt 1
t=linspace(-pi,pi,65)
t=t[:-1]
dt=t[1]-t[0]
fmax=1/dt
y=sin(sqrt(2)*t)
y[0]=0 #the centrepoint is set to 0
y=fftshift(y) #make y start from y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65)
w=w[:-1]

subplot(2,1,1)
plot(w, abs(Y), lw=2)
xlim([-10,10])
ylabel(r"$|Y| ->$", size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)$")
grid()
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$ ->",size=16)
xlabel(r"$\omega$ ->",size=16)
grid() 
show()

#plotting the function over multiple time periods
t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
plot(t, sin(sqrt(2)*t), color='red')
plot(t2,sin(sqrt(2)*t2), color='blue')
plot(t3,sin(sqrt(2)*t3), color='blue')
xlabel('time ->')
ylabel(r'$\sin\left(\sqrt{2}t\right)$ ->')
title(r'plot of $\sin(\sqrt{2}t)$')
grid()
show()

#replicating the blue portion alone over all the periods
plot(t,sin(sqrt(2)*t),'ro')
plot(t2,sin(sqrt(2)*t),'bo')
plot(t3,sin(sqrt(2)*t),'bo')
xlabel('time ->')
ylabel('y ->')
title(r'plot of $\sin(\sqrt{2}t)$ wrapping every 2$\pi$')
grid() 
show()


#DFT of a ramp
z=t
z[0]=0
z=fftshift(z)
Z=fftshift(fft(z))/64
semilogx(abs(w),20*log10(abs(Z)),lw=2)
xlim([1,10])
ylim([-20,0])
xticks([1,2,5,10],["1","2","5","10"],size=16) #cool feature puts markers on the x-axis
ylabel(r"$|Y|$ (dB)",size=16)
title(r"Spectrum of a digital ramp")
xlabel(r"$\omega$",size=16) 
grid() 
show()

#testing the hamming window
n=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*n/63))
y=sin(sqrt(2)*t1)*wnd
plot(t1,y,'bo',lw=2)
plot(t2,y,'ro',lw=2)
plot(t3,y,'ro',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)\times w(t)$ with $t$ wrapping every $2\pi$ ") 
grid()
show()

#DFT of this sequence
y[0]=0 # the sample corresponding to -tmax should be set zero
y=fftshift(y) # make y start with y(t=0) 
Y=fftshift(fft(y))/64.0 
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-8,8]) 
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2) 
xlim([-8,8])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

#using more samples to improve accuracy
t3=linspace(-4*pi,4*pi,257);t3=t3[:-1]
dt1=t3[1]-t3[0]
fmax1=1/dt1
n=arange(256)
wnd=fftshift(0.54+0.46*cos(2*pi*n/256))
y=sin(sqrt(2)*t3) 
y=y*wnd 
y[0]=0 # the sample corresponding to -tmax should be set zero
y=fftshift(y) # make y start with y(t=0) 
Y=fftshift(fft(y))/256.0 
w=linspace(-pi*fmax1,pi*fmax1,257);w=w[:-1]
subplot(2,1,1)
plot(w,abs(Y),'b',w,abs(Y),'bo',lw=2) 
xlim([-4,4])
ylabel(r"$|Y|$",size=16) 
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$") 
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16) 
grid(True) 
show()

#worked examples end here
def selector(t,n,w=None,d=None):     #selector function 
    if(n==1):
        return pow(cos(0.86*t),3)
    elif(n==2):
        return cos(16*(1.5+t/(2*pi))*t)
    elif(n==3):
        if(w!=None and d!=None):
            return cos(w*t+d)
    else:
            return pow(cos(0.86*t),3)

def window_fn(n,N):
    return (0.54+0.46*cos(2*pi*n/N))


'''
Function to find Discrete Fourier Transform
Arguments:
 low_lim,up_lim -> lower & upper limit for time
 no_points      -> Sampling rate
 n              -> function number
 wo             -> frequency
 d              -> phase difference
 noise_const    -> Gaussian contribution
 '''

def findFFT(low_lim,up_lim,no_points,n,window=True,wo=None,d=None,noise_const=0):
    t = linspace(low_lim,up_lim,no_points+1)[:-1]
    dt=t[1]-t[0]
    fmax=1/dt
    N = no_points
    y = selector(t,n,wo,d)+noise_const*randn(len(t)) #noise_const is only given if required.Randn() emulates noise
    
    if(window):
        n1=arange(N)
        wnd=fftshift(window_fn(n1,N))
        y=y*wnd
        
    y[0]=0        # the sample corresponding to -tmax should be set zero
    y=fftshift(y) # make y start with y(t=0)
    Y=fftshift(fft(y))/N

    w = linspace(-pi*fmax,pi*fmax,N+1)[:-1]        
    return t,Y,w
'''
Function to plot Magnitude and Phase spectrum for given function
Arguments:
 t              -> time vector
 Y              -> DFT computed
 w              -> frequency vector
 Xlims,Ylims    -> limits for x&y axis for spectrum
 plot_title -> title of plot
'''

def plot_FFT(t,Y,w,Xlims,plot_title,dotted=False,Ylims=None):
    
    subplot(2,1,1)
    
    if(dotted):
        plot(w,abs(Y),'bo',lw=2)
    else:
        plot(w,abs(Y),'b',lw=2)

    xlim(Xlims)
        
    ylabel(r"$|Y(\omega)| \to$")
    title(plot_title)
    grid(True)
    
    subplot(2,1,2)
    ii=where(abs(Y)>0.005)
    plot(w[ii],angle(Y[ii]),'go',lw=2)

    if(Ylims!=None):
        ylim(Ylims)
    
    xlim(Xlims)
    ylabel(r"$\angle Y(j\omega) \to$")
    xlabel(r"$\omega \to$")
    grid(True)
    show()
    
#Question 2 for cos^3((Wo)t)
t,Y,w = findFFT(-4*pi,4*pi,256,1,False)
Xlims = [-8,8]
plot_title = r"Spectrum of $\cos^3(\omega_o t)$ without windowing"
plot_FFT(t,Y,w,Xlims,plot_title)

t,Y,w = findFFT(-4*pi,4*pi,256,1,True)
Xlims = [-8,8]
plot_title = r"Spectrum of Windowed $\cos^3(\omega_o t)$"
plot_FFT(t,Y,w,Xlims,plot_title)

#Question 3 and 4 estimating w and delta
def estimator(noise_const=0):   #function estimates the values of w and delta
    w_actual = np.random.uniform(0.5,1.5) #generates a random omega value
    delta_actual = (randn())              #random delta value
    t,Y,w = findFFT(-1*pi,1*pi,128,3,True,w_actual,delta_actual,noise_const)
    
    ii = where(w>=0)
    w_cal = sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2)   #do weighted average over positive w
    i = abs(w-w_cal).argmin() #find w at which it is closest to the estimated w0
    delta_cal = angle(Y[i])   #find delta at that w0 since it'll be the phase
    print("Calculated value of w0 (to 4 decimal places) ",round(w_cal,4))
    print("Calculated value of delta (to 4 decimal places) ",round(delta_cal,4))
    print("Actual value of w0 (to 4 decimal places) ",round(w_actual,4))
    print("Actual value of delta (to 4 decimal places)",round(delta_actual,4))

print("without noise")
estimator() #without noise
print("with noise")
estimator(0.1) #with noise

#section 5 and 6
def chirp(t):
    return cos(16*(1.5+t/(2*pi))*t)
t = linspace(-pi,pi,1025)[:-1]
dt=t[1]-t[0]
fmax=1/dt
N = 1024
y = chirp(t)
Y=fftshift(fft(y))/N
w = linspace(-pi*fmax,pi*fmax,N+1)[:-1]       
Xlims = [-100,100]
plot_title = r"Spectrum of Non-Windowed Chirped Signal"
plot_FFT(t,Y,w,Xlims,plot_title)
t = linspace(-pi,pi,1025)[:-1]
dt=t[1]-t[0]
fmax=1/dt
N = 1024
y = chirp(t)

n1=arange(N)
wnd=fftshift(window_fn(n1,N))
y=y*wnd
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/1024.0
w = linspace(-pi*fmax,pi*fmax,1025)[:-1]       

Xlims = [-100,100]
plot_title = r"Spectrum of Windowed Chirped Signal"
plot_FFT(t,Y,w,Xlims,plot_title)

#the code for the time frequency surface plot(6th)
t_array = split(t,16)  #splitting time into 16
Y_mag = zeros((16,64))
Y_phase = zeros((16,64))

for i in range(len(t_array)):
    n = arange(64)
    wnd = fftshift(0.54+0.46*cos(2*pi*n/64))
    y = cos(16*t_array[i]*(1.5 + t_array[i]/(2*pi)))*wnd
    y[0]=0
    y = fftshift(y)
    Y = fftshift(fft(y))/64.0
    Y_mag[i] = abs(Y)
    Y_phase[i] = angle(Y)

t = t[::64]	
w = linspace(-fmax*pi,fmax*pi,64+1); w = w[:-1]
t,w = meshgrid(t,w)

fig1=figure(1)
ax = fig1.add_subplot(111, projection='3d')
surf=ax.plot_surface(w,t,Y_mag.T,cmap='viridis',linewidth=0, antialiased=False)
fig1.colorbar(surf, shrink=0.5, aspect=5)
ax.set_title('surface plot');
ylabel(r"$\omega\rightarrow$")
xlabel(r"$t\rightarrow$")
show()
