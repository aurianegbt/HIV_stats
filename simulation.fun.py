import numpy as np
import matplotlib.pyplot as plt

class Simulation(object):
    def __init__(self,T_max,p,y0):
        ''' Initialise une simulation pour les paramètres T_max, y0 et p '''
        self.T = T_max
        self.p = p
        self.ic = y0

        (S0,I0,R0)=y0
        self.N = sum(y0)

        self.S = [y0[0]/sum(y0)]
        self.I = [y0[1]/sum(y0)]
        self.R = [y0[2]/sum(y0)]
        self.listt=[0]
        self.det_times=[]

    def taux(self):
        ''' taux de saut du modèle SIR pour S,I,R fixé '''
        (lambda_0,lambda_1,lambda_2,lambda_3,mu_0,mu_1,c) = self.p
        S = self.S[-1]
        I = self.I[-1]
        R = self.R[-1]

        N=self.N
        return(np.array([lambda_0,
        mu_0*S*N,
        lambda_1*S*I*N,
        mu_1*I*N,
        lambda_2*I*N,
        lambda_3*I*S*N]
        ))


    def run(self):
        ''' Génère une simulation du modèle SIR sur [0;T_max], pour les paramètres p et
        condition initiale y0=(S0,I0,R0) '''
        (lambda_0,lambda_1,lambda_2,lambda_3,mu_0,mu_1,c) = self.p
        t_current = self.listt[-1]
        (S0,I0,R0)=self.ic #nbr d'individu pas taux ici

        SIR_current=np.array([S0,I0,R0])


        while t_current < self.T and SIR_current[1]>0:
            # On calcule la date du prochain évènement :
            e = np.random.exponential(scale=sum(taux(self)))
            t_next = t_current + e

            # On détermine le type de saut :
            u = np.random.uniform(low=0,high=sum(self))
            cum_taux = np.cumsum(self)
            R_aux=self.R[-1]
            #On met à jours la liste des temps
            self.listt.append(t_next)
            self.S.append(self.S[-1])
            self.I.append(self.I[-1])
            self.R.append(self.R[-1]*np.exp(-c*e))

            self.listt.append(t_next)
            # On effectue le saut :
            if u< cum_taux[0] :
                SIR_current += np.array([1,0,0])
                self.N = sum(SIR_current)
                self.S.append(self.S[-1] +1/self.N)
                self.I.append(self.I[-1])
                self.R.append(self.R[-1])
            elif cum_taux[0]<=u and u<cum_taux[1]:
                SIR_current += np.array([-1,0,0])
                self.N = sum(SIR_current)
                self.S.append(self.S[-1] -1/self.N)
                self.I.append(self.I[-1])
                self.R.append(self.R[-1])
            elif cum_taux[1]<=u and u<cum_taux[2]:
                SIR_current += np.array([-1,+1,0])
                self.S.append(self.S[-1] - 1/self.N)
                self.I.append(self.I[-1] + 1/self.N)
                self.R.append(self.R[-1])
            elif cum_taux[2]<=u and u<cum_taux[3]:
                SIR_current += np.array([0,-1,0])
                self.N = sum(SIR_current)
                self.S.append(S[-1])
                self.I.append(I[-1]-1/self.N)
                self.R.append(R[-1])
            elif cum_taux[3]<=u and u<cum_taux[4]:
                SIR_current += np.array([0,-1,+1])
                self.S.append(S[-1])
                self.I.append(I[-1]-1/self.N)
                self.R.append(R[-1]+1/self.N)
                self.det_times.append([t_next,1]) #1 pour détection spontanée
            elif cum_taux[4]<=u and u<cum_taux[5]:
                SIR_current += np.array([0,-1,+1])
                self.S.append(self.S[-1])
                self.I.append(self.I[-1]-1/self.N)
                self.R.append(self.R[-1]+1/self.N)
                self.det_times.append([t_next,2]) #2 pour détection par contact tracing
            else:
                self.S.append(self.S[-1])
                self.I.append(self.I[-1])
                self.R.append(self.R[-1])
            t_current = t_next
        return(self)

    def simuR(self):
        ''' Fonction pour extraire R^1_t, R^2_t de la simulation sim '''
        times_det = self.det_times[:,0]
        type_det = self.det_times[:,1]
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

    def plotR1R2(self):
        (listt,R1,R2)=simu_data(self)

        plt.plot(listt,R1,label="$R^1(t)$",color='gold')
        plt.plot(listt,R2,label="$R^2(t)$",color='firebrick')

        plt.legend()
        plt.title('Simulation du système SIR, $R^1, R^2$')

        plt.show()

    def plotsimu(self):
        ''' Plot et renvoie une smiluation du modèle SIR sur [0;T_max], pour les
        paramètres p et condition initiale y0=(S0,I0,R0) '''


        plt.plot(self.listt,self.S,label="$s(t)$",color='gold')
        plt.plot(self.listt,self.I,label="$i(t)$",color='firebrick')
        plt.plot(self.listt,self.R,label="$r(t)$",color='limegreen')

        plt.legend()
        plt.title('Simulation du système SIR')
