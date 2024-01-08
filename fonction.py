import numpy as np

###FONCTIONS NUMERIQUES###
def RK3(f, Y0, temps, alpha, theta):
    c1 = 1 - 1/(2*alpha) - (3*alpha-2)/(6*alpha*theta)
    c2 = 1/(2*alpha) + (3*alpha-2)/(6*alpha*(theta-alpha))
    c3 = -(3*alpha-2)/(6*theta*(theta-alpha))
    a32 = -theta*(theta-alpha)/(alpha*(3*alpha-2))
    a31 = theta - a32
    n = len(Y0)
    p = len(temps)
    Y = np.zeros((n, p))
    Y[:,0] = Y0
    for i in range(p-1):
        h = temps[i+1] - temps[i]
        s1 = temps[i]
        s2 = temps[i] + alpha*h
        s3 = temps[i] + theta*h
        k1 = Y[:,i]
        k2 = Y[:,i] + h*alpha*f(s1, k1)
        k3 = Y[:,i] + h*(a31*f(s1,k1) + a32*f(s2,k2))
        Y[:,i+1] = Y[:,i] + h*(c1*f(s1,k1) + c2*f(s2,k2) + c3*f(s3,k3))
    return Y

def e(X,Y,N):
    E=[]
    for i in range(0,N):
        E.append(np.linalg.norm(X[:,i]-Y[:,i]))
    return(max(E))

def f_osc(t,Y):
    x,y = Y
    return np.array([-y,x])

def f_raide(t,Y):
    x,y=Y
    return np.array([-100*x,-y])

def f_bernoulli(t,Y):
    y=Y[0]
    return(np.array([(t+1)*np.sqrt(y)-2*y]))

def f_riccati(t,Y):
    return np.array([(-Y[0]**3)/2])

###FONCTIONS D'AFFICHAGE###

#Cette fonction va permettre d'afficher le menu dans le terminal
def afficher_menu():
    print("Menu:")
    print("1. Résoudre à l'aide de RK3 l'équation de l'oscillateur harmonique")
    print("2. Résoudre à l'aide de RK3 une équation de Bernoulli")
    print("3. Résoudre à l'aide de RK3 une équation raide")
    print("4. Résoudre à l'aide de RK3 une équation de Riccati")
    print("5. Récapitulatif des erreurs en fonction du nb de pts")
    print("6. Quitter")

#Cette fonction va permettre de rentrer les paramètres au clavier
def acquisition(b):
    alpha=1
    theta=1
    #On fait l'acquisition de paramètres différents selon b (b=0 pour les choix 1 à 4, b=1 pour le tracé du tableau des erreurs)
    if b==0:
        N=int(input("choisissez le nomre de points N:"))
        t0=float(input("choisissez t0:"))
        tf=float(input("choisissez tf:"))
        alpha,theta=choix_param()
        return(N,t0,tf,alpha,theta)

    elif b==1:
        t0=float(input("choisissez t0:"))
        tf=float(input("choisissez tf:"))
        alpha,theta=choix_param()
        return(t0,tf,alpha,theta)
    
#On choisit de faire les calculs suivant un schéma proposé ou en sélectionnant soit-même les paramètres
def choix_param():
    print("1: Choix arbitraire des paramètres")
    print("2: Méthode de Ralston")
    print("3: Méthode de Heun")
    print("4: Méthode de Kutta")
    print("5: Méthode de SSPRK3")
    choix=int(input("Choisissez la méthode de résolution: "))
    alpha=1
    theta=1
    if choix==1:
        alpha_temp=input("choisissez alpha:")
        theta_temp=input("choisissez theta:")
        if '/' in alpha_temp and '/' not in theta_temp:
            num,den=alpha_temp.split('/')
            num=float(num)
            den=float(den)
            alpha=num/den
            theta=float(theta_temp)
        elif '/' in theta_temp and '/' not in alpha_temp:
            num,den=theta_temp.split('/')
            num=float(num)
            den=float(den)
            theta=num/den
            alpha=float(alpha_temp)
        elif '/' in alpha_temp and '/' in theta_temp:
            num,den=alpha_temp.split('/')
            num=float(num)
            den=float(den)
            alpha=num/den
            num,den=theta_temp.split('/')
            num=float(num)
            den=float(den)
            theta=num/den     
        else:
            alpha=float(alpha_temp)
            theta=float(theta_temp)
    if choix==2:
        alpha=1/2
        theta=3/4
    if choix==3:
        alpha=1/3
        theta=2/3
    if choix==4:
        alpha=1/2
        theta=1
    if choix==5:
        alpha=1
        theta=1/2
    return(alpha,theta)