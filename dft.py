from math import pi, sqrt, degrees
from math import sin,cos,atan2
import cmath


def amplitudeResponse(N,A0,k,m):
    """
    Calculate the amplitude response of DFT with given input
    
    :param N: Number of points of input
    :param A0: Peak value of input sinusoid
    :param k: Number of cycles of input sinusoid
    :param m: The DFT output index
    :returns: Returns the amplitude at specified index
    """
    sinc = 0
    v = (k - m)
    x = pi * v
    if v == 0:
        sinc = 1.0
    else:
        sinc = ( sin(x) / (x))
        
    return ((A0 * N)/2) * sinc  



def amplitudeResponseRange(k, N=32, A0=1.0):
    amp = []
    for m in range(0,N):
        a = amplitudeResponse(N, A0, k, m)
        amp.append( a )
                    
    return amp


def dft(x):
    result = []
    N = len(x);
    pi2 = 2 * pi
    for m in range(0,N):
        sumr = 0.0
        sumi = 0.0
        for n in range(0,N):
            sumr = sumr + (x[n] * cos( (pi2 * n * m) / N))
            sumi = sumi + (x[n] * sin( (pi2 * n * m) / N))
        result.append(complex(round(sumr,3)+0, round(sumi,3)+0))
    
    return result

def magnitude(z):
    return round(sqrt( pow(z.real, 2) + pow(z.imag,2)))

def phase(z):
    return round(degrees(atan2(z.imag, z.real)))


def windowTriangle(x):
    n = 0
    N = len(x)
    pivotf = N/2.0
    pivoti = int(round(pivotf))
    
    while n < pivoti:
        coef = n/pivotf
        x[n] = x[n] * coef
        n = n + 1

    n = pivoti
    while n < N:
        coef = 2 - (n/pivotf)
        x[n] = x[n] * coef
        n = n + 1

def windowHanning(x):
    n = 0
    N = len(x)
    pi2 = 2 * pi
    while n < N:
       x[n] = x[n] * (0.5 - (0.5 * cos( pi2 * n / N) ))
       n = n + 1
                            
def windowHamming(x):
    n = 0
    N = len(x)
    pi2 = 2 * pi
    while n < N:
       x[n] = x[n] * (0.54 - (0.46 * cos( pi2 * n / N)  )) 
       n = n + 1
       
       
def timeForSpectralRes(resF):
    return 1.0/resF
    