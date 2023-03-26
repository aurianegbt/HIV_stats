import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.style.use('ggplot')
from simulation.fun import Simulation

## Fonctions utiles
def Kd(t,delta):
    ''' Noyau D'Epanechinov sur [-delta;delta] '''
    a = 4/3
    if np.abs(t)<=delta:
        return(a/delta*(1-(t/delta)**2))
    else:
        return(0)

def Norme_fun(sim,data):
    ''' Calcul la norme L1 de la simulation et des données '''
    (times_det,R1,R2)=simu_data(sim)
    (times_det_obs,R1_obs,R2_obs)=data

    if max(times_det_obs)<=max(times_det):
        R1 = interp1d(times_det,R1)(times_det_obs)
        R2 = interp1d(times_det,R2)(times_det_obs)
        n=len(times_det_obs)
    else:
        R1_obs = interp1d(times_det_obs,R1_obs)(times_det)
        R2_obs = interp1d(times_det_obs,R2_obs)(times_det)
        n=len(times_det)


    N1 = [np.abs(R1[i]-R1_obs[i]) for i in range(n)]
    N2 = [np.abs(R2[i]-R2_obs[i]) for i in range(n)]
    return((np.mean(N1),np.mean(N2)))

def calibration_delta(listsimu,data):
    ''' Calibration du delta dans Kd pour ne garder que les 1% plus proches simulations'''
    # On a fixé P_delta=0.01, là on a une liste listN de nombre entre 0 et +oo et on veut garder les 1% plus petites, delta sera le quantile
    listN1=[]
    listN2=[]
    N=len(listsimu)
    for sim in listsimu:
        (N1,N2)=Norme_fun(sim,data)
        listN1.append(N1)
        listN2.append(N2)
    ind1 = np.argsort(listN1)
    delta1 = max(listsimu[ind1[:(N//100)]])
    ind2 = np.argsort(listN2)
    delta2 = max(listsimu[ind2[:(N//100)]])
    return(delta1,delta2)


def sample_prior(self):
    lambda_0 = 0
    mu_0 = 0
    mu_1 = 0

    lambda_1 = np.random.uniform(low=0.5,high=0.7)
    lambda_2 = np.random.uniform(low=0.01,high=0.3)
    lambda_3 = np.random.uniform(low=0.01,high=0.3)
    c = np.random.uniform(low=0.01,high=0.3)

    p=(lambda_0,lambda_1,lambda_2,lambda_3,mu_0,mu_1,c)

    self.p = p
    return(p)

## Simulation des données
T_max= 1000000
y0=(1000,10,0)

lambda_0 = 0
mu_0 = 0
mu_1 = 0
lambda_1 = 0.6
lambda_2 = 0.2
lambda_3 = 0.2
c = 0.12         # Taux pour le contact tracing
p=(lambda_0,lambda_1,lambda_2,lambda_3,mu_0,mu_1,c)

Sim = Simulation(T_max,p,y0)
Sim.run()

data = Sim.simuR()

## Simulations

N=10000
listsimu=[]
listprior=[]
for i in range(N):
    p_prior = sample_prior()
    listprior+=[p_prior]
    listsimu+=[Simulation(T_max,p_prior,y0)]

delta1,delta2=calibration_delta(listsimu,data)

for sim in listsimu:
    print(Norme_fun(sim,data))

listW=[]
for sim in listsimukeep:
    (N1,N2) = Norme_fun(sim,data)
    listW.append(Kd(N1,delta)*Kd(N2,delta))

posterior_l1=[]
posterior_l2=[]
posterior_l3=[]
posterior_c=[]
for i in range(N):
    (lambda_0i,lambda_1i,lambda_2i,lambda_3i,mu_0i,mu_1i,ci)=listprior[i]
    posterior_l1+=[lambda_1i]
    posterior_l2+=[lambda_2i]
    posterior_l3+=[lambda_3i]
    posterior_c+=[ci]


## Tracer de la posterior

# Lambda_1
plt.hist(posterior_l1,weights=listW,bins=30)
plt.show()

# Lambda_2
plt.hist(posterior_l2,weights=listW,bins=30)
plt.show()

# Lambda_3
plt.hist(posterior_l3,weights=listW,bins=30)
plt.show()

# c
plt.hist(posterior_c,weights=listW,bins=30)
plt.show()