import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

def taux(S,I,R,N,p):
    ''' taux de saut du modèle SIR pour S,I,R fixé '''
    (lambda_0,lambda_1,lambda_2,lambda_3,mu_0,mu_1,c) = p
    return(np.array([
        lambda_0,
        mu_0*S*N,
        lambda_1*S*I*N,
        mu_1*I*N,
        lambda_2*I*N,
        lambda_3*I*S*N
    ]))

def simu(T_max,p,y0):
    ''' Génère une simulation du modèle SIR sur [0;T_max], pour les paramètres p et
    condition initiale y0=(S0,I0,R0) '''
    (lambda_0,lambda_1,lambda_2,lambda_3,mu_0,mu_1,c) = p
    t_current = 0
    listt=[t_current]
    det_times=[]
    (S0,I0,R0)=y0
    
    SIR_current=np.array([S0,I0,R0])
    N = sum(SIR_current)
    S = [SIR_current[0]/N]
    I = [SIR_current[1]/N]
    R = [SIR_current[2]/N]
    
    
    while t_current < T_max and SIR_current[1]>0:
        # On calcule la date du prochain évènement :
        e = np.random.exponential(scale=sum(taux(S[-1],I[-1],R[-1],N,p)))
        t_next = t_current + e
        
        # On détermine le type de saut : 
        u = np.random.uniform(low=0,high=sum(taux(S[-1],I[-1],R[-1],N,p)))
        cum_taux = np.cumsum(taux(S[-1],I[-1],R[-1],N,p))
        R_aux=R[-1]
        #On met à jours la liste des temps
        listt.append(t_next)
        S.append(S[-1])
        I.append(I[-1])
        R.append(R[-1]*np.exp(-c*e))
        
        listt.append(t_next)
        # On effectue le saut : 
        if u< cum_taux[0] :
            SIR_current += np.array([1,0,0])
            N = sum(SIR_current)      
            S.append(S[-1] +1/N)
            I.append(I[-1])
            R.append(R[-1])
        elif cum_taux[0]<=u and u<cum_taux[1]:
            SIR_current += np.array([-1,0,0])
            N = sum(SIR_current)      
            S.append(S[-1] -1/N)
            I.append(I[-1])
            R.append(R[-1])
        elif cum_taux[1]<=u and u<cum_taux[2]:
            SIR_current += np.array([-1,+1,0])   
            S.append(S[-1] - 1/N)
            I.append(I[-1] + 1/N)
            R.append(R[-1])
        elif cum_taux[2]<=u and u<cum_taux[3]:
            SIR_current += np.array([0,-1,0]) 
            N = sum(SIR_current)  
            S.append(S[-1])
            I.append(I[-1]-1/N)
            R.append(R[-1])
        elif cum_taux[3]<=u and u<cum_taux[4]:
            SIR_current += np.array([0,-1,+1]) 
            S.append(S[-1])
            I.append(I[-1]-1/N)
            R.append(R[-1]+1/N)
            det_times.append([t_next,1]) #1 pour détection spontanée
        elif cum_taux[4]<=u and u<cum_taux[5]:
            SIR_current += np.array([0,-1,+1]) 
            S.append(S[-1])
            I.append(I[-1]-1/N)
            R.append(R[-1]+1/N)
            det_times.append([t_next,2]) #2 pour détection par contact tracing
        else:
            S.append(S[-1])
            I.append(I[-1])
            R.append(R[-1])
        t_current = t_next
    return(listt,S,I,R,det_times)

def plot_simu(T_max,p,y0):
    ''' Plot et renvoie une smiluation du modèle SIR sur [0;T_max], pour les
     paramètres p et condition initiale y0=(S0,I0,R0) '''
    (listt,S,I,R,det_times)= simu(T_max,p,y0)


    plt.plot(listt,S,label="$s(t)$",color='gold')
    plt.plot(listt,I,label="$i(t)$",color='firebrick')
    plt.plot(listt,R,label="$r(t)$",color='limegreen')

    plt.legend()
    plt.title('Simulation du système SIR')

    plt.show()

    return(listt,S,I,R,det_times)

def plot_simu2(simu):
    ''' Plot et renvoie une smiluation du modèle SIR sur [0;T_max], pour les
     paramètres p et condition initiale y0=(S0,I0,R0) '''
    (listt,S,I,R,det_times)= simu


    plt.plot(listt,S,label="$s(t)$",color='gold')
    plt.plot(listt,I,label="$i(t)$",color='firebrick')
    plt.plot(listt,R,label="$r(t)$",color='limegreen')

    plt.legend()
    plt.title('Simulation du système SIR')

    plt.show()
    
def simu_data(sim):
    ''' Fonction pour extraire R^1_t, R^2_t de la simulation sim '''
    det_times=np.array(sim[4])
    times_det = det_times[:,0]
    type_det = det_times[:,1]
    N=len(type_det)
    listtR=[0]
    R1=[0]
    R2=[0]
    for n in range(N):
        listtR.append(times_det[n])
        if type_det[n]==1:
            R1.append(R1[-1]+1)
            R2.append(R2[-1])
        else :
            R1.append(R1[-1])
            R2.append(R2[-1]+1)
    return(listtR,R1,R2)

def plot_R1R2(sim):
    (listt,R1,R2)=simu_data(sim)
    
    plt.plot(listt,R1,label="$R^1(t)$",color='gold')
    plt.plot(listt,R2,label="$R^2(t)$",color='firebrick')

    plt.legend()
    plt.title('Simulation du système SIR, $R^1, R^2$')

    plt.show()
    
    
def sample_prior():
    ''' Tirage des paramètres $\theta$ selon les priors'''
    # Lambda_0,mu_0 est supposé nulle pour les simulations 
    lambda_0 = 0
    mu_0 = 0
    mu_1 = 0
    
    lambda_1 = np.random.uniform(low=0.5,high=0.7)
    lambda_2 = np.random.uniform(low=0.01,high=0.3)
    lambda_3 = np.random.uniform(low=0.01,high=0.3)
    c = np.random.uniform(low=0.01,high=0.3)
    
    p=(lambda_0,lambda_1,lambda_2,lambda_3,mu_0,mu_1,c)
    return(p)