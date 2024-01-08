import numpy as np
import matplotlib.pyplot as plt
import fonction as fct
import timeit

while True:
    fct.afficher_menu()
    choix = int(input("Entrez votre choix (1-6): "))

    if choix==1:
            N,t0,tf,alpha,theta=fct.acquisition(0)
            
            #Calcul de la solution numérique
            Y0 = np.array([1, 0])
            LesTemps = np.linspace(t0, tf, N)

            t_debut = timeit.default_timer() # Mesure du temps de calcul
            LesY_RK3_osc = fct.RK3(fct.f_osc, Y0, LesTemps,alpha,theta)
            t_fin = timeit.default_timer() #Calcul du temps écoulé
            duree = t_fin-t_debut
            plt.plot(LesY_RK3_osc[0, :], LesY_RK3_osc[1, :],color='red',linestyle="dashed", label="Approximation numérique")
            
            #Calcul de la solution exacte
            Sol_exact=np.array([np.cos(LesTemps),np.sin(LesTemps)])
            plt.plot(Sol_exact[0,:],Sol_exact[1,:],color='purple',alpha=0.5, label="Solution exacte")
            
            #Affichage
            plt.axis('equal')
            plt.title("Comparaison solution exacte/approximation numérique")
            plt.xlabel("y_1(t)")
            plt.ylabel("y_2(t)")
            plt.legend()
            plt.show()

            #Calcul de l'erreur et du coût de calcul en temps
            erreur=fct.e(Sol_exact,LesY_RK3_osc,N)
            print(f"l'erreur est {erreur}")
            print(f"Le temps de calcul est de {duree}")


    if choix==2:
            N,t0,tf,alpha,theta=fct.acquisition(0)

            #Calcul de la solution numérique
            LesTemps = np.linspace(t0, tf, N)
            Y0=np.array([1])

            t_debut = timeit.default_timer() # Mesure du temps de calcul
            LesY_RK3_bernoulli=fct.RK3(fct.f_bernoulli,Y0,LesTemps,alpha,theta)
            t_fin = timeit.default_timer() #Calcul du temps écoulé
            duree = t_fin-t_debut
            plt.plot(LesTemps,LesY_RK3_bernoulli[0],color='red',linestyle="dashed", label="Approximation numérique")
            
            #Calcul de la solution exacte
            Sol_exact=(LesTemps/2 +np.exp(-LesTemps))**2
            Sol_exact = Sol_exact.reshape(1,N)
            plt.plot(LesTemps,Sol_exact[0],color='orange',alpha=0.5, label="Solution exacte")

            #Affichage
            plt.legend()
            plt.axis('equal')
            plt.title("Comparaison solution exacte/approximation numérique")
            plt.xlabel("t")
            plt.ylabel("y(t)")
            plt.show()
            
            #Calcul de l'erreur et du coût de calcul en temps
            erreur=fct.e(LesY_RK3_bernoulli,Sol_exact,N)
            print(f"l'erreur est {erreur}")
            print(f"Le temps de calcul est de {duree}")


    if choix==3:
            N,t0,tf,alpha,theta=fct.acquisition(0)

            #Calcul de la solution numérique
            Y0=np.array([1,1])
            LesTemps = np.linspace(t0, tf, N)

            t_debut = timeit.default_timer() # Mesure du temps de calcul
            LesY_RK3_raide=fct.RK3(fct.f_raide,Y0,LesTemps,alpha,theta)
            t_fin = timeit.default_timer() #Calcul du temps écoulé
            duree = t_fin-t_debut
            plt.plot(LesY_RK3_raide[0],LesY_RK3_raide[1],color="red",linestyle="dashed", label="Approximation numérique")
            
            # Calcul de la solution exacte
            Sol_exact=np.array([np.exp(-100*LesTemps),np.exp(-LesTemps)])
            plt.plot(Sol_exact[0],Sol_exact[1],color="green",alpha=0.5, label="Solution exacte")
            
            #Affichage
            plt.title("Comparaison solution exacte/approximation numérique")
            plt.legend()
            plt.xlabel("y_1(t)")
            plt.ylabel("y_2(t)")
            plt.show()

            #Calcul de l'erreur et du coût de calcul en temps
            erreur=fct.e(Sol_exact,LesY_RK3_raide,N)
            print(f"l'erreur est {erreur}")
            print(f"Le temps de calcul est de {duree}")

    if choix==4:
            N,t0,tf,alpha,theta=fct.acquisition(0)

            #Calcul de la solution numérique
            Y0=np.array([1])
            LesTemps=np.linspace(t0,tf,N)

            t_debut = timeit.default_timer() # Mesure du temps de calcul
            LesY_RK3_riccati=fct.RK3(fct.f_riccati,Y0,LesTemps,alpha,theta)
            t_fin = timeit.default_timer() #Calcul du temps écoulé
            duree = t_fin-t_debut
            plt.plot(LesTemps,LesY_RK3_riccati[0],color="red",linestyle="dashed",label="Approximation numérique")  

            Sol_exact=np.array([1/np.sqrt(1+LesTemps)])
            Sol_exact=Sol_exact.reshape(1,N)
            plt.plot(LesTemps,Sol_exact[0],label="Solution exacte",alpha=0.5)

            #Affichage
            plt.legend()
            plt.title("Comparaison solution exacte/approximation numérique")
            plt.xlabel("t")
            plt.ylabel("y(t)")
            plt.show()

            #Calcul de l'erreur et du coût de calcul en temps
            erreur=fct.e(Sol_exact,LesY_RK3_riccati,N)
            print(f"l'erreur est {erreur}")
            print(f"Le temps de calcul est de {duree}")
    
    if choix==5:
        while True:
                print("Choisissez l'équation différentielle:")
                print("1. L'équation de l'oscillateur harmonique")
                print("2. Une équation de Bernoulli")
                print("3. Une équation raide")
                print("4. Une équation de Riccati")
                print("5. Quitter")
                choix_recap = int(input("Entrez votre choix (1-5): "))


                LesN=[101,201,401,501,801,1001,2001,4001,8001,10001]
                nb_pts=len(LesN)
                Tab_Erreur=np.zeros((nb_pts,2))
                Tab_Erreur[:,0]=LesN

                if choix_recap==1:
                    t0,tf,alpha,theta=fct.acquisition(1)
                    for i in range(nb_pts):
                        Y0 = np.array([1, 0])
                        LesTemps_recap = np.linspace(t0, tf, LesN[i])
                        LesY_RK3_osc_recap = fct.RK3(fct.f_osc, Y0, LesTemps_recap,alpha,theta)
                        Sol_exact_recap=np.array([np.cos(LesTemps_recap),np.sin(LesTemps_recap)])
                        Tab_Erreur[i,1]=fct.e(Sol_exact_recap,LesY_RK3_osc_recap,LesN[i])

                    print("Pour l'équation de l'oscillateur harmonique, on a")
                    print(Tab_Erreur)

                elif choix_recap==2:
                    t0,tf,alpha,theta=fct.acquisition(1)
                    for i in range(nb_pts):
                        Y0 = np.array([1])
                        LesTemps_recap = np.linspace(t0, tf, LesN[i])
                        LesY_RK3_bernoulli_recap = fct.RK3(fct.f_bernoulli, Y0, LesTemps_recap,alpha,theta)
                        Sol_exact_recap=(LesTemps_recap/2 +np.exp(-LesTemps_recap))**2
                        Sol_exact_recap=Sol_exact_recap.reshape(1, LesN[i])
                        Tab_Erreur[i,1]=fct.e(Sol_exact_recap,LesY_RK3_bernoulli_recap,LesN[i])

                    print("Pour l'équation de Bernoulli, on a")
                    print(Tab_Erreur)                       

                elif choix_recap==3:
                    t0,tf,alpha,theta=fct.acquisition(1)         
                    for i in range(nb_pts):
                        Y0 = np.array([1, 1])
                        LesTemps_recap = np.linspace(t0, tf, LesN[i])
                        LesY_RK3_raide_recap = fct.RK3(fct.f_raide, Y0, LesTemps_recap,alpha,theta)
                        Sol_exact_recap=np.array([np.exp(-100*LesTemps_recap),np.exp(-LesTemps_recap)])
                        Tab_Erreur[i,1]=fct.e(Sol_exact_recap,LesY_RK3_raide_recap,LesN[i])

                    print("Pour l'équation raide, on a")
                    print(Tab_Erreur)  
                
                elif choix_recap==4:
                    t0,tf,alpha,theta=fct.acquisition(1)
                    for i in range(nb_pts):
                        Y0 = np.array([1])
                        LesTemps_recap = np.linspace(t0, tf, LesN[i])
                        LesY_RK3_riccati_recap = fct.RK3(fct.f_riccati, Y0, LesTemps_recap,alpha,theta)
                        Sol_exact_recap=np.array([1/np.sqrt(1+LesTemps_recap)])
                        Tab_Erreur[i,1]=fct.e(Sol_exact_recap,LesY_RK3_riccati_recap,LesN[i])

                    print("Pour l'équation de Riccati, on a")
                    print(Tab_Erreur) 

                elif choix_recap==5:
                        break



    if choix==6:
            break  # Sortir de la boucle lorsque l'option "Quitter" est choisie





