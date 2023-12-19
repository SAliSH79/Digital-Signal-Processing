#!/usr/bin/env python
# coding: utf-8

# # DSP_Lab Project
# # 9723042
# # 3 - 1
# SeyedAli SeyedHosseini

# In[259]:


import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
# Import libraries


# ### 3 - 1 ) a

# In[260]:


# Initialization
R = 0.8
f0 = 500
fs = 10**4
# test
print(fs)


# In[261]:


w0 = (2*np.pi)*f0/fs # Discrete Frequency
print(w0)


# In[262]:


a1 = -2*R*np.cos(w0) #a1
print(a1)


# In[263]:


a2 = R**2 #a2
print(a2)


# In[264]:


G = (1 - R)*np.sqrt((1 - 2*R*np.cos(2*w0) + R**2)) #Coeff of Numerator
print(G)


# In[265]:


# Frequency Response
num = [G]
den = [1 ,a1 ,a2]
w ,h = signal.freqz(num, den, worN=1024, whole = True) #fs = 10**4) #Frequency Response


# In[266]:


f = w/np.pi #Normalization
fig, ax1 = plt.subplots()
ax1.plot(f,np.abs(h)**2,'-r')
ax1.grid()
ax1.set_xlabel('freq (normilized by pi)')
ax1.set_title("Frequency Response |H(ejw)|")
ax1.set_ylabel("Mag")
ax1.legend("|H(ejw)|" , loc = 'upper center')
plt.xlim([0, 2])
plt.ylim([-0.25 , 1.25])
plt.show()


# ### 3 - 1 ) b

# In[267]:


n = np.arange(300)
#print(n)
h1 = G/np.sin(w0)
#print(h1)
h2 = R ** n
#print(h2)
h3 = np.sin((n + 1)*w0)
#print(h3)
h = h1*h2*h3 #impulse Response in Formula

N = 300
x = np.zeros(N)
x[0] = 1 #Creating Impulse in 2nd Sample
#print(x)
#print("\n")
hf = signal.lfilter(num,den,x)
#print(hf)

y = np.zeros(N) #preallocation
for i in range (1 , N) :
    if i == 1 :
        y[i] = G*1
    elif i == 2 :
        y[i] = -a1 * y[i - 1]
    else : 
        y[i] = -a1*y[i - 1] -a2*y[i - 2]
#print(y)


# In[268]:


## PLOTTING 

plt.figure
plt.plot(h,"*r",hf,"--k",y,"b")
#plt.stem(y,"b") isn`t working
plt.legend(('syntax','Convolution','Recursive'),loc='best')
plt.grid()
plt.title("Impulse Response");
plt.xlabel("time");
plt.ylabel("Amp");
plt.xlim([-5, 305]);
plt.ylim([-0.04 , 0.22]);


# ### 3 - 1 ) c

# In[269]:


s = np.cos(w0*n) #input signal
#print(s)
#print("\n")
v = np.random.randn(1,N) #random noise Generation
#print(v)
x = s + v #Input of Filter
#print("\n")
#print(X)
y = np.zeros(N)
w1 = 0
w2 = 0
for i in range(1,N) :
    y[i] = -a1*w1 - a2*w2 + G*x[0,i] # x is 2d dimension because of np.cos of n which is a vector too
    w2 = w1
    w1 = y[i]


# In[270]:


## PLOTTING 

plt.figure
plt.plot(y,"g",s,"r")
plt.legend(('FIlter Output',"S[n]"),loc='best')
plt.grid()
plt.title("Output of H(z) - S[n]");
plt.xlabel("time");
plt.ylabel("Amp");
plt.xlim([-5, 305]);
plt.ylim([-2.2 ,1.7]);


# In[271]:


#print(n)


# In[272]:


#print(s)


# ### 3 - 1 ) d

# In[273]:


n = np.arange(300)
v = np.random.randn(1,N) #random noise Generation
y_v = np.zeros(N)
w1 = 0
w2 = 0
for i in range(1,N) :
    y_v[i] = -a1*w1 - a2*w2 + G*v[0,i] # V is 2d dimension because of np.cos of n which is a vector too
    w2 = w1
    w1 = y[i]


# In[274]:


## PLOTTING 

plt.figure
plt.plot(y_v,"k",v[0],"--r")
plt.legend(('FIlter Noise Output',"V[n]"),loc='best')
plt.grid()
plt.title("Noisy Output of H(z) - V[n]");
plt.xlabel("time");
plt.ylabel("Amp");
#plt.xlim([-5, 305]);
#plt.ylim([-2.2 ,1.7]);


# In[275]:


#Plotting

fig,(ax1,ax2) = plt.subplots(1,2)
fig.suptitle('Noise Plots before and after Filter ');
ax1.plot(v[0],"--r");
ax2.plot(y_v,"k");
ax1.grid();
ax1.set_xlabel('time(n)');
ax1.set_title("V[n]");
ax1.set_ylabel("Amp");
ax1.legend("V[n]" , loc = 'upper center');
ax2.grid();
ax2.set_xlabel('time(n)');
ax2.set_title("y_v[n]");
ax2.set_ylabel("Amp");
ax2.legend("y_v[n]" , loc = 'best');


# ### 3 - 1)e

# In[276]:


s1v = np.var(v[0])
print("Variance of Noise is",s1v)
s1yv = np.var(y_v)
print("Variance of y_v (Output noise from filter) is",s1yv)
NRR = s1yv/s1v  #Noise Reduction Ratio is var(y_v)/var(v)
print('So NRR is', NRR) 


# In[277]:


# Now lets check Formula
Num = (1 + 2*R**2) #Correct Formula
Den = (1 + R)*(1 + 2*R*np.cos(2*w0) + R**2) #Correct Formula
NRR_frml = Num/Den
print('So NRR by formula is', NRR_frml)


# ## The End
