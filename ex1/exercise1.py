import math

def simpson(f, a, b, n):
    """Approximate the Integral of f in the interval
    between a and b using n equidistant sample points using 
    Simpsons rule. The number of sample points can be even or odd."""
    # Enforce oddity
    n = max (n, 2)
    h = (b-a)/(n-1)
    odd_n = n -1 + n%2 # Enforce oddity
    S = f(a) + f(a+(odd_n-1)*h)
    # Integrate over alternating coefficients
    for i in range(1, n-1, 2):
        S += 4*f(a + i * h)
    for i in range(2, n-2, 2):
        S += 2*f(a + i * h)
    S *= h/3
    # Case n even
    if (n%2 == 0):
        S += (h/12)* (5*f(b) + 8*f(b-h) - f(b-2*h))
    return S

# Lets output a test value
simpson(math.sin, 0, math.pi/2, 8)

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# Theoretical error bound
def deltaSimpson(a, b, n, k):
    return k*((b-a)**5)/(180*(n**4))

# Collect Simpson integrals and errors for different n in [0, pi/2]
n = []
S = [] 
ErrBound = []  # theoretical error
# Interation Number, Number of Sample Points, S_T(h_n), Delta S 
print("{:10}| {:10}| {:20}| {:20}".format("Iteration", "N", "S_s", "Delta S"))
for k in range(1, 8):
    n.append(10**k)
    S.append(simpson(math.sin, 0, math.pi/2, n[k-1]))
    ErrBound.append(deltaSimpson(0, math.pi/2, n[k-1], 1))
    print("{:10}| {:10}| {:20}| {:20}".format(k, n[k-1], S[k-1], ErrBound[k-1]))

# Actual error is 1 - S_T(h)
Deviation = [] # actual error
for i in S:
    print("1 - {} = {}".format(i, abs(1-i)))
    Deviation.append(abs(1-i))

line1, = plt.plot(n, Deviation, '--', label='Actual Error')
line2, = plt.plot(n, ErrBound, '-.', label='Theoretical Error')
# Log-Log plot of actual error versus theoretical error
plt.grid()
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Number of samples N")
plt.legend(handles=[line1, line2], loc='best')
plt.show()

def romberg(f, a, b, n):
    """Romberg integration of f in [a,b] with n levels."""
    h = b-a # Interval size   
    Table = np.zeros((n, n))   
    Table[0][0] = ((h/2) * (f(a) + f(b))) # Insert S_T(h_0)
    # Trapezoidal rule
    for k in range(1,n):
        h /= 2
        # Collect new sample points
        new_samples = 0
        for i in range(1, 2**k, 2):
            new_samples += f(a+i*h)
        Table[k][0] = (Table[k-1][0]/2) + h*new_samples
        
    # Richardson extrapolation
    for c in range(1, n):
        pow4 = (4**c) - 1
        for r in range(c, n):
            Table[r][c] = Table[r][c-1] + ((Table[r][c-1] - Table[r-1][c-1])/pow4) 
          
    return Table[n-1][n-1]

# Test value, 10 steps should get us around 20 digits of accurary
romberg(math.sin, 0, math.pi/2, 10)

k = 22
# Tabulate the R_ii for the three test functions for increasing values of i
Rii = np.zeros((3, k-2))

# go get a coffee
for i in range(2, k):
    Rii[0][i-2] = romberg((lambda x: math.exp(x)),      0, 1,         i)
    Rii[1][i-2] = romberg((lambda x: math.sin(8*x)**4), 0, 2*math.pi, i)
    Rii[2][i-2] = romberg((lambda x: math.sqrt(x)),     0, 1,         i)
    
# calculate deviation from analytical result
Delta = np.zeros((3, k-2))
for i in range(0, k-2):
    Delta[0][i] = abs((math.exp(1)-1) - abs(Rii[0][i]))
    Delta[1][i] = abs(3*math.pi/4 - abs(Rii[1][i]))
    Delta[2][i] = abs((2/3) - abs(Rii[2][i]))

x = np.linspace(0,k-2,k-2)

fig, ((ax11, ax12, ax13), (ax21, ax22, ax23)) = plt.subplots(2, 3)
fig.set_size_inches((12, 6))
fig.suptitle('Romberg integration and the deviations from the analytical result for increasing stepnumber i')

ax11.plot(x,Rii[0:1].T,'r--')
ax11.set_title("Rii Equation 1")
ax12.plot(x,Rii[1:2].T,'g--')
ax12.set_title("Rii Equation 2")
ax13.plot(x,Rii[2:3].T,'b--')
ax13.set_title("Rii Equation 3")
ax21.plot(x,Delta[0:1].T,'r--')
ax21.set_title("Delta Eq1")
ax22.plot(x,Delta[1:2].T,'g--')
ax22.set_title("Delta Eq2")
ax23.plot(x,Delta[2:3].T,'b--')
ax23.set_title("Delta Eq3")
plt.setp((ax21, ax22, ax23), yscale="log")
plt.show()
