import pandas as pd
import numpy as np
import math
from tqdm import tqdm
import matplotlib.pyplot as plt

# Load the CSV file
data = pd.read_csv('/Users/shaswataroy/Documents/Jing/Xenium/JAK-STAT Pathway/yamada/stat1n_dimer_data.csv')

# Access the data
time = data['time_hours']
stat1n_dimer = data['STAT1n_star_STAT1n_star']

# Time Steps
dt = .1
tspan = np.arange(0,T,dt)

# Parameters
k15 =  0.019
k16 =  0.025
k17 =  0.022
tau_avg_2 = 22
q2 = 4./tau_avg_2


# Kernel of Gamma function
def Gam(tau, p, q):
    return ((q**p)/math.gamma(p))*tau**(p-1)*np.exp(-q*tau)

def response(t,q,n):
    integral_response = 0
    t_i = round(t/dt)

    for i in range(t_i):
        tau = i*dt
        integral_response += Gam(tau,4,q)*stat1n_dimer[t_i-i]*dt
    
    return integral_response

def model(x,t):
    
    i = int(round(t/dt))
    
    v = k15 + k16*response(t,q2,4) - k17*stat1n_dimer[i]

    return v

# Time step
tspan  = np.linspace(0, 100, 1000)

irf1 = np.zeros(len(tspan))
irf1[0] = 0

for i in tqdm(range(len(tspan)-1)):
    irf1[i+1] = model(irf1[i],tspan[i])*dt + irf1[i]

plt.plot(tspan, stat1n_dimer)
plt.xlabel('Time (min)')
plt.ylabel('STAT1n-Dimer')
plt.title('STAT1n-Dimer vs Time')
plt.show()
