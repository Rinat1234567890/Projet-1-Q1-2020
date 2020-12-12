"""
Source: Pecheur Charles 2018.
"""

import math as m
import matplotlib.pyplot as plt
import numpy as np

"""
***Variables de base: la barge est supposÃ©e stable, carrÃ©e et n'oscille pas.***
"""

long   = 0.60                  #Longueur de la barge. (m)
haut   = 0.115                  #Hauteur de la plateforme. (m)
dist_0 = 0.00                  #Distance initiale de la charge. (m)
dist_t = 0.40                  #Distance finale de la charge. (m)
dist_v = 0.01                  #Vitesse a laquelle la barge se deplace. (m/s)
coeff_amort=0.6                #Coefficient d'amortissement

mass_tot=6.3                   #Masse totale de la grue et plateforme (kg) 
mass_obj=0.130                 #Masse de l'objet deplace (kg)

### Centre de gravite. Il faut les calculer Ã  la main.

cg_x_0 = 0.10
cg_y_0 = 0.15

###### Simulation

### Constantes

g = 9.81         # gravitation [m/s**2]
rho_eau = 1000   # masse volumique [kg/m**3]

"""
*****Parametres de la simulation*****
"""

dt = 0.001           # pas (dt) [s]
end = 30             # duree [s]
theta_0 = 0.0        # position initiale [rad]
omega_0 = 0.0        # vitesse initiale [rad/s]

t = np.arange(0, end, dt)
theta = np.empty_like(t)          #position angulaire [rad]
omega = np.empty_like(t)          #vitesse angulaire [rad/s]
alpha = np.empty_like(t)          #acceleration initiale [rad/s**2]
dist = np.empty_like(t)

soulev= np.empty_like(t)
submer= np.empty_like(t)
soulev2= np.empty_like(t)
submer2= np.empty_like(t)
haut_enf=mass_tot/(rho_eau*long**2)

### Fonction de rotation. Elle met Ã  jour les valeurs nÃ©cessaires pour le calcul.

def rotation(i):
    global delta_h
    global cg_x_final
    global cg_y_final
    global cp_x
    global cp_y
    global cp_x_final
    global cp_y_final
    
    delta_h= (long/2) * m.tan(theta[i])
    
    cg_x=cg_x_0+(dist[i]*mass_obj)/mass_tot
    cg_y=cg_y_0                                         #definition du centre de poussee
    
    cg_x_final=cg_x*m.cos(theta[i])-cg_y*m.sin(theta[i]) 
    cg_y_final=cg_x*m.sin(theta[i])+cg_y*m.cos(theta[i]) #rotation du centre de gravite

    cp_x = (delta_h*long) / (6*haut_enf)
    cp_y = haut_enf/2 - (delta_h**2)/(6*haut_enf)       #definition du centre de poussee
    
    cp_x_final=cp_x*m.cos(theta[i])-cp_y*m.sin(theta[i])
    cp_y_final=cp_x*m.sin(theta[i])+cp_y*m.cos(theta[i]) #rotation du centre de poussee

if haut_enf>=haut:
    print("Attention: haut_enf > haut. La barge est submergee")

def simulation():
    """
    post: calcule l'oscillation de la barge
    """
    theta[0] = theta_0
    omega[0] = omega_0
    dist[0]  = dist_0
    
    global delta_h
    global cg_x_final
    global cg_y_final
    global cp_x
    global cp_y
    global cp_x_final
    global cp_y_final    

    for i in range(len(t)-1):
        rotation(i)
        mf=(cg_x_final-cp_x_final)*mass_tot*g
        mf_amort=-omega[i]*g*coeff_amort
        mf_tot=mf+mf_amort
        
        alpha[i]=mf_tot/mass_tot #Ceci calcule le couple de redressement
        omega[i+1]=omega[i]+alpha[i]*dt
        theta[i+1]=theta[i]+omega[i]*dt
        alpha[i+1]=alpha[i] 
        #Ces trois conditions correspondent au calcul par iteration de l'angle, pulsation et acceleration angulaires dans le modele physique.
        if abs(dist[i]-dist_t)>=abs(2*dist_v*dt): #Cette condition 'if' calcule le deplacement de la charge en fonction du temps.
            dist[i+1]=dist[i]+dist_v*dt
        else:
            dist[i+1]=dist_t 
    haut_max=m.atan(2*(haut-haut_enf)/long)
    delta_h_max=m.atan(2*(haut_enf)/long)
    #Calcul des angles correspondant aux situations critiques du bateau (soulevement et submersion). Leurs noms sont un peu trompeurs.
    
    for i in range(len(t)):
        soulev[i]=(delta_h_max)
        submer[i]=(haut_max)
        soulev2[i]=-(delta_h_max)
        submer2[i]=-(haut_max)
        #Dessin des points critiques
        
def graphiques(): 
    """
    Cette fonction dessine et legende les graphiques.
    """
    plt.figure(1) #Creation d'une image blanche sur laquelle se trouveront les graphiques
    plt.subplot(3,1,1) #Creation du graphique indiquant l'evolution de l'angle au cours du temps
    plt.plot(t,theta, label="angle")
    plt.plot(t,soulev,color='orange',linestyle='dashed')
    plt.plot(t,submer,color='red',linestyle='dashed') 
    #Les trois .plot() dessinent les angles du bateau et les angles maximales dans le premier graphique
    
    print("Orange == Soulevement\nRouge == Submersion")
    plt.legend() #Commande qui legende le graphique.
    plt.subplot(3,1,2) #Creation du deuxieme graphique: celui ci indique l'evolution de la pulsation angulaire au cours du temps
    plt.plot(t,omega, label="v. angulaire")
    plt.legend()
    plt.subplot(3,1,3) #creation du troisieme graphique: celui ci indique l'evolution de l'acceleration angulaire au cours du temps
    plt.plot(t,alpha, label="a. angulaire")
    plt.legend()
    plt.show() #Renvoie l'image avec les graphiques

simulation()
graphiques()
#Etant donne que ce fichier n'est pas cense etre importe, la condition 'if __name__ == "__main__"' n'a pas ete introduite.

print("Angle final:",theta[29999]*180/m.pi,"degres")
print("Soulevement:",soulev[29999]*180/m.pi,"degres")
print("Submersion:",submer[29999]*180/m.pi,"degres") 
#Les trois print() renvoient les dernieres valeurs pour l'angle theta et les angles maximales.
