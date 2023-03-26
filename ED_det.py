import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

def euler(f,y0,listt,p=None):

    """ Renvoie la solution de l'EDO y'=f(p,t,y) sur [0;T] par la méthode d'Euler,
     de condition initiale (t0,y0), avec p les paramètres du modèle, et listt la
     subdivision de [0;T] utilisée   """

    m=len(y0)

    sol = np.zeros((len(listt),m),dtype=object)
    sol[0]=y0

    for n in range(0,len(listt)-1):
        tn  = listt[n]
        dtn = listt[n+1]-tn
        yn  = sol[n]

        sol[n+1]=yn+dtn*f(tn,yn,p)

    return(sol)

def f(t,y,p):
    ''' Second membre du modèle SIR considérée'''
    (lbd0,lbd1,lbd2,lbd3,mu0,mu1,c)=p

    (s,i,r)=y

    fy = np.array([ lbd0 - mu0*s - lbd1*s*i,
          lbd1*s*i - (mu1+lbd2)*i - lbd3*i*r,
          lbd2*i + lbd3*i*r - c*r ])

    return(fy)

def sol_ed(listt,p,y0):
    ''' Renvoie la solution du système limite du modèle SIR '''
    return(euler(f,y0,listt,p))

def plot_sol_ed(listt,p,y0):
    ''' Plot de la solution du système limite du modèle SIR '''
    plt.figure(figsize=(10.,5.))

    sol=sol_ed(listt,p,y0)

    plt.plot(listt,sol[:,0],label="$s(t)$",color='gold')
    plt.plot(listt,sol[:,1],label="$i(t)$",color='firebrick')
    plt.plot(listt,sol[:,2],label="$r(t)$",color='limegreen')

    plt.legend()
plt.title('Solution du système ($S$) par la méthode d\'Euler,\n avec $(\lambda_0,\lambda_1,\lambda_2,\lambda_3,\mu_0,\mu_1,,c)$='+str(p))

    plt.show()