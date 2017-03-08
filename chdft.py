from dft import dft, sin,pi, magnitude as mag, phase, windowTriangle,\
    windowHanning, windowHamming,timeForSpectralRes
from dft import amplitudeResponseRange as ampRange
import matplotlib.pyplot as plt

def sectionDftWavePoint(t):
    return sin( 2 * pi * 1000 * t) + 0.5 *(sin( 2 * pi * 2000*t+(3*pi/4)))

def waveSine3Hz(t):
    return sin( 2 * pi * 3000 * t)
    
    
def sectionDftWave(f, N = 8, fs = 8000):
    delta = 1.0/fs
    
    result = []
    for t in range(0,N):
        result.append( f( t * delta) )
    
    return result

def sectionLeakage():
    # DFT Leakage
    amp8Hz = ampRange(8.00)     # No leakage
    amp85Hz = ampRange(8.50)    # Leakage
    amp875Hz = ampRange(8.75)   # Leakage
        
    plt.figure().suptitle("8.00Hz")
    posamp8Hz = list(filter(lambda a: a >= 0, amp8Hz))
    plt.bar(range(0,14), posamp8Hz[:14],width=1/32,label="8Hz", align="center")
    plt.xlabel("DFT index")
    plt.ylabel("Amplitude Response")
    plt.grid(True)
    
    plt.figure().suptitle("8.5Hz")
    posamp875Hz = list(filter(lambda a: a >= 0, amp875Hz))
    plt.bar(range(0,14), posamp875Hz[:14],width=1/32,label="8.75Hz", align="center", color="red")
    plt.xlabel("DFT index")
    plt.ylabel("Amplitude Response")
    plt.grid(True)

    plt.figure().suptitle("8.75Hz")
    posamp = list(filter(lambda a: a >= 0, amp875Hz))
    plt.bar(range(0,14), posamp[:14],width=1/32,label="8.75Hz", align="center")
    plt.xlabel("DFT index")
    plt.ylabel("Amplitude Response")
    plt.grid(True)
    

    
    plt.show()
    
def sectionDftInteger():
    wave = sectionDftWave(sectionDftWavePoint, 8, 8000)
    result = dft(wave)
    
    x = 0
    for z in result:

        za = z 
        # Get the correct sign for polar form
        # magnitude, phase
        if z.imag != 0:
            za = complex(z.real, z.imag * -1)
        else:
            za = z

        # Change sign for displaying resolved rectagular form
        # X(m) = a - jb
        op = "-"
        if z.imag < 0.0:
            op = "+"
            z = complex(z.real, z.imag * -1)
            
            
        
        print("x({})\t{:.3f} {} j{:.3f}\t\t{}, {}".format(x,z.real, op, z.imag, mag(za), phase(za)))
        x = x + 1
        
    mags = list(map(lambda z: mag(z), result))
    fig, ax = plt.subplots()
    
    
    ax.bar(range(0,len(mags)), mags, align="center", color="green", ecolor="black", width=1/8)
    ax.set_xticks(range(0,len(mags)))
    plt.show()
    
def sectionDftNonInteger():
    wave = sectionDftWave(sectionDftWavePoint, 8, 8500)
    result = dft(wave)
    
    x = 0
    for z in result:

        za = z 
        # Get the correct sign for polar form
        # magnitude, phase
        if z.imag != 0:
            za = complex(z.real, z.imag * -1)
        else:
            za = z

        # Change sign for displaying resolved rectagular form
        # X(m) = a - jb
        op = "-"
        if z.imag < 0.0:
            op = "+"
            z = complex(z.real, z.imag * -1)
            
            
        
        print("x({})\t{:.3f} {} j{:.3f}\t\t{}, {}".format(x,z.real, op, z.imag, mag(za), phase(za)))
        x = x + 1
        
    mags = list(map(lambda z: mag(z), result))
    fig, ax = plt.subplots()
    
    
    ax.bar(range(0,len(mags)), mags, align="center", color="green", ecolor="black", width=1/8)
    ax.set_xticks(range(0,len(mags)))
    plt.show()
    
def sectionDft():
    sectionDftNonInteger() 
    
    
def sectionWindows():
    size = 100
    x = [1.0] * size 
    windowTriangle(x)
    
    plt.figure().suptitle("Triangle window response")
    plt.bar(range(0,len(x)), x,width=1/32, align="center")
    plt.xlabel("W(n)")
    plt.ylabel("Coefficient Value")
    plt.grid(True)
    
    x = [1.0] * size
    windowHanning(x)
    plt.figure().suptitle("Hanning window response")
    plt.bar(range(0,len(x)), x,width=1/32, align="center")
    plt.xlabel("W(n)")
    plt.ylabel("Coefficient Value")
    plt.grid(True)
    
    x = [1.0] * size
    windowHamming(x)
    plt.figure().suptitle("Hamming window response")
    plt.bar(range(0,len(x)), x,width=1/32, align="center")
    plt.xlabel("W(n)")
    plt.ylabel("Coefficient Value")
    plt.grid(True)
    
    plt.show()
        
def zeroPad(x, N):
    for i in range(0,N):
        x.append(0.0)    
def sectionZeroPadding():
    N = 32
    wave = sectionDftWave(waveSine3Hz, N, 8000)
    result = dft(wave)
    

    mags = list(map(lambda z: mag(z), result))
    plt.figure().suptitle("L=16 N=16")
    plt.bar(range(0,len(mags)), mags, align="center", color="green", ecolor="black", width=1/8)
   
    wave = sectionDftWave(waveSine3Hz, N, 8000)
    zeroPad(wave, N * 1)
    result = dft(wave)
    mags = list(map(lambda z: mag(z), result))
    plt.figure().suptitle("L=16 N=32")
    plt.bar(range(0,len(mags)), mags, align="center", color="green", ecolor="black", width=1/8)
    
    wave = sectionDftWave(waveSine3Hz, N, 8000)
    zeroPad(wave, N * 2)
    result = dft(wave)
    mags = list(map(lambda z: mag(z), result))
    plt.figure().suptitle("L=16 N=64")
    plt.bar(range(0,len(mags)), mags, align="center", color="green", ecolor="black", width=1/8)

    
    wave = sectionDftWave(waveSine3Hz, N, 8000)
    zeroPad(wave, N * 4)
    result = dft(wave)
    mags = list(map(lambda z: mag(z), result))
    plt.figure().suptitle("L=16 N=128")
    plt.bar(range(0,len(mags)), mags, align="center", color="green", ecolor="black", width=1/8)
    
    
    wave = sectionDftWave(waveSine3Hz, N, 8000)
    zeroPad(wave, N * 8)
    result = dft(wave)
    mags = list(map(lambda z: mag(z), result))
    plt.figure().suptitle("L=16 N=256")
    plt.bar(range(0,len(mags)), mags, align="center", color="green", ecolor="black", width=1/8)
    plt.show()
    
def sectionSpectralResolutionSize():
    print("For 3000Hz resolution: {} seconds".format(timeForSpectralRes(3000)))
    print("For 300Hz resolution: {} seconds".format(timeForSpectralRes(300)))
    print("For 30Hz resolution: {} seconds".format(timeForSpectralRes(30)))
    print("For 3Hz resolution: {} seconds".format(timeForSpectralRes(3)))