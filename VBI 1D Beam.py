# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 20:20:11 2026
VBI Analysis with a bridge modelled as a 1D beam 

@author: Brillant
"""
import numpy as np
from scipy.linalg import eigh, lu_factor, lu_solve

gamma, beta =  0.5, 0.25     # Methode de Newmark à acceleration constante (inconditionnelement stable)
#nbr_elements = 3            # Nombre de subdivision de chaque travée
p = mw = Mv = kv = cv = kB = v = a0 = a1 = a2 = a3 = a4 = a5 = a6 = a7 = 0
epsilon_1 = epsilon_2 = 0   # Coefficient d'amortissement modal du materiau (beton par ex)
g = 9.81   #m/s^2
# Cette fonction retourne la matrice de rigité elementaire de chaque element        
def Matrice_rig_el_pout_1D(E,I,L):
    return (E*I/L**3)*np.array([[12, 6*L, -12, 6*L],
                                [6*L, 4*L**2, -6*L, 2*L**2],
                                [-12, -6*L, 12, -6*L],
                                [6*L, 2*L**2, -6*L, 4*L**2]])

# Cette fonction retourne la matrice de masse elementaire de chaque element        
def Matrice_masse_el_pout_1D(rho,A,L):
    return (rho*A*L/420)*np.array([[156, 22*L, 54, -13*L],
                                [22*L, 4*L**2, 13*L, -3*L**2],
                                [54, 13*L, 156, -22*L],
                                [-13*L, -3*L**2, -22*L, 4*L**2]])

# Cette fonction defini les fonctions de formes(fonction d'interpolation d'hermite)
def N(x, l):
    N1 = 1-3*(x/l)**2+2*(x/l)**3
    N2 = 3*(x/l)**2-2*(x/l)**3
    N3 = x-2*x**2/l+x**3/l**2
    N4 = -x**2/l+x**3/l**2
    return np.array([[N1],[N3],[N2],[N4]])

# Cette matrice recherche l'element sur lequel la force ponctuelle est appliquée
# position_force est un vecteur ligne (ex: np.array([x1,x2,x3]))
def Localise_element(position_force, Tab_connect):
    numero_element = np.zeros([np.shape(position_force)[0]], dtype = int)
    position_locale_de_la_force = np.zeros([np.shape(position_force)[0]])
    for j in range(np.shape(position_force)[0]):
        elt = 0
        continuer = True
        while(continuer):
            if position_force[j] < 0 or position_force[j] > Tab_connect[-1, 6] :
                numero_element[j] = -1  #-1 signifie que la force n'est pas sur la poutre
                position_locale_de_la_force[j] = -1
                continuer = False
            elif position_force[j] >= Tab_connect[elt, 5] and position_force[j] <= Tab_connect[elt, 6] : 
                numero_element[j] = elt
                position_locale_de_la_force[j] = position_force[j] - Tab_connect[elt, 5]
                continuer = False
            else:
                elt += 1
                if elt >= len(Tab_connect):
                    numero_element[j] = -1
                    position_locale_de_la_force[j] = -1
                    break
    return numero_element, position_locale_de_la_force     # Retourne -1 si la force n'est pas sur la poutre                   

# Cette fonction permet à l'utilisateur de definir les donnéés d'entrer
def Donnee_entree():
    mw = float(input("\n Masse de la roue du vehicule (en kg): mw = "))
    Mv = float(input(" Masse de la caisse du vehicule (en kg): Mv = "))
    kv = float(input(" Raideur de suspension (en N/m): kv = "))
    cv = float(input(" Amortissement de suspension (en N.s/m): cv = "))
    kB = float(input(" Raideur de contact (Ballast ou Pneu) (en N/m): kB = "))
    v = float(input(" Vitesse du vehicule (en m/s): V = "))
    epsilon_1 = float(input(" Coefficient d'amortissement modal du materiau: xi 1 = "))
    epsilon_2 = float(input(" Coefficient d'amortissement modal du materiau: xi 2 = "))
    nbr_elements = int(input(" Vous souhaitez mailler chaque travé en combien d'élément: Nombre d'élément = "))
    nbr_dt = int(input(" Nombre de Pas de temps: Nombre de Delta t = "))
    tol = float(input(" Tolérance des resultats du calcul (exemple: 10^-3) : "))
    max_it = int(input(" Nombre maximal d'iteration (exemple: 30) : "))
    print("\n ---------- La durée de l'etude sera defini de sorte à permettre au vehicule de traverser entierement le pont ----------")
    return mw, Mv, kv, cv, kB, v, nbr_elements, nbr_dt, epsilon_1, epsilon_2, tol, max_it

# Cette fonction definie les irregularité de la route
def r(x):
    irregularite = 0*x
    return irregularite



# Cette fonction permet à l'utilisateur de definir les condition aux limites et renvoi les ddl a supprimé    
def Numeros_DDL_a_supprimer(nombre_total_de_noeuds):
    ddl_supprimer = []
    print("\n Le nombre de travée determine le nombre d'appuis intermediaire. \n Chosissez un numero correspondant aux conditions aux limites des extremités gauche et droite de la poutre (qui ne doit pas etre hypostatique) : ")
    n = int(input(" 1. Appuis-Appuis \n 2. Appuis-Libre \n 3. Libre-Appuis \n 4. Encastré-Encastré \n 5. Encastré-Appuis \n 6. Appuis-Encastré \n 7. Encastré-Libre \n 8. Libre-Encastré \n 9. Libre-Libre (à eviter dans le cas du pont) \n\n Choix N°: "))
    if n == 1:
        ddl_supprimer.append(0)
        ddl_supprimer.append(nombre_total_de_noeuds*2-2)
    elif n == 2:
        ddl_supprimer.append(0)
    elif n == 3:
        ddl_supprimer.append(nombre_total_de_noeuds*2-2)
    elif n == 4:
        ddl_supprimer.append(0)
        ddl_supprimer.append(1)
        ddl_supprimer.append(nombre_total_de_noeuds*2-1)
        ddl_supprimer.append(nombre_total_de_noeuds*2-2)
    elif n == 5:
        ddl_supprimer.append(0)
        ddl_supprimer.append(1)
        ddl_supprimer.append(nombre_total_de_noeuds*2-2)
    elif n == 6:
        ddl_supprimer.append(0)
        ddl_supprimer.append(nombre_total_de_noeuds*2-1)
        ddl_supprimer.append(nombre_total_de_noeuds*2-2)
    elif n == 7:
        ddl_supprimer.append(0)
        ddl_supprimer.append(1)
    elif n == 8:
        ddl_supprimer.append(nombre_total_de_noeuds*2-1)
        ddl_supprimer.append(nombre_total_de_noeuds*2-2)
    elif n == 9:
        ddl_supprimer = []
    else:
        print("Votre choix ne peux pas etre pris en compte. La suite des resultats n'ai donc pas adapté à votre cas.")
    return ddl_supprimer  # Numeros des DDL à supprimer aux extremité

# Cette fonction defini la geometrie de la poutre et retourne sa table de connectivité 
def Definition_pout(nbr_element):
    position_element = 0
    nbr_travee = int(input(" Entrez le nombre de travée: "))
    longueur_travées = np.zeros([nbr_travee])
    inertie_travées =  np.zeros([nbr_travee])
    module_young_travées = np.zeros([nbr_travee])
    masse_volumique_travées =  np.zeros([nbr_travee])
    section_travées =  np.zeros([nbr_travee])
    Tab_connect = np.zeros([nbr_travee*nbr_element, 9])
    for i in range(nbr_travee):
        exposant = ("ère")*(i==0)+("ème")*(i!=0)
        longueur_travées[i] = float(input(" Entrer la portée de la "+str(i+1)+" "+str(exposant)+" travée (en m): L"+str(i+1)+" = "))
        inertie_travées[i] = float(input(" Entrer le moment d'inertie de la "+str(i+1)+" "+str(exposant)+" travée (en m^4): I"+str(i+1)+" = "))
        module_young_travées[i] = float(input(" Entrer le module de Young de la "+str(i+1)+" "+str(exposant)+" travée (en MPa): E"+str(i+1)+" = "))*10**6
        masse_volumique_travées[i] = float(input(" Entrer la masse volumique de la "+str(i+1)+" "+str(exposant)+" travée (en kg/m^3): Rho"+str(i+1)+" = "))
        section_travées[i] = float(input(" Entrer la section de la "+str(i+1)+" "+str(exposant)+" travée (en m^2): A"+str(i+1)+" = "))
        subdivision = np.linspace(0, longueur_travées[i], nbr_element+1)    # nbre de point = nbre element + 1
        l = subdivision[1]-subdivision[0]
        # tab_connect contient repectivement pour chaque element: 0)Noeud de depard, 1)Noeud d'arriver, 2)longueur, 3)inetie, 4)module de young,  5)position de debut, 6)position de fin, 7)masse volumique de la travée, 8)section de la travée 
        for elm in range(nbr_element):
            j = elm + nbr_element*i
            Tab_connect[j, 0] = j
            Tab_connect[j, 1] = j+1
            Tab_connect[j, 2] = l
            Tab_connect[j, 3] = inertie_travées[i]
            Tab_connect[j, 4] = module_young_travées[i]
            Tab_connect[j, 5] = position_element + l*elm # debut
            Tab_connect[j, 6] = position_element + l*(elm+1) # fin
            Tab_connect[j, 7] = masse_volumique_travées[i]
            Tab_connect[j, 8] = section_travées[i]
        position_element += longueur_travées[0*(i==0)+(i)*(i!=0)] 
    #Matrice de rigidité élémentaire des éléments  discrétisés 
    Mat_elem=[] # liste devrant stocker les matrices de rig. élementaire
    Masse_elm = [] # liste devrant stocker les matrices de masse élementaire
    for elt in range(np.shape(Tab_connect)[0]):
        k_elem = Matrice_rig_el_pout_1D(Tab_connect[elt,4],
                                         Tab_connect[elt,3],
                                         Tab_connect[elt,2])
        Mat_elem.append(k_elem)
        m_elem = Matrice_masse_el_pout_1D(Tab_connect[elt,7],
                                         Tab_connect[elt,8],
                                         Tab_connect[elt,2])
        Masse_elm.append(m_elem)
    ddl_intermediaire_supprimé = []    
    for i in range(nbr_travee-1):
        nbr_ddl_travee_precedente = ((nbr_element+1)*2-2)*i
        ddl_intermediaire_supprimé.append(nbr_ddl_travee_precedente + (nbr_element+1)*2-2)
    
    nombre_total_de_noeuds = nbr_travee*(nbr_element+1)-(nbr_travee - 1)
    return Tab_connect, Mat_elem, Masse_elm, nombre_total_de_noeuds, ddl_intermediaire_supprimé

# Cette fonction fait l'assemblage des matrices de rigidité ou de masse
# ddl_supprimer = ddl_intermediaire_supprimé + ddl_supprimer(aux extremités) (concatenation des deux listes)
def Assemblage_matrice(nombre_total_de_noeuds, Tab_connect, Mat_elem, ddl_supprimer):
    k_global = np.zeros([nombre_total_de_noeuds*2,nombre_total_de_noeuds*2])
    for i in range(np.shape(Tab_connect)[0]):
        k_global[2*(int(Tab_connect[i,0])):2*(int(Tab_connect[i,0]))+2,2*(int(Tab_connect[i,0])):2*(int(Tab_connect[i,0]))+2]=k_global[2*(int(Tab_connect[i,0])):2*(int(Tab_connect[i,0]))+2,2*(int(Tab_connect[i,0])):2*(int(Tab_connect[i,0]))+2]+Mat_elem[i][0:2,0:2]
        k_global[2*(int(Tab_connect[i,1])):2*(int(Tab_connect[i,1]))+2,2*(int(Tab_connect[i,1])):2*(int(Tab_connect[i,1]))+2]=k_global[2*(int(Tab_connect[i,1])):2*(int(Tab_connect[i,1]))+2,2*(int(Tab_connect[i,1])):2*(int(Tab_connect[i,1]))+2]+Mat_elem[i][2:4,2:4]
        k_global[2*(int(Tab_connect[i,0])):2*(int(Tab_connect[i,0]))+2,2*(int(Tab_connect[i,1])):2*(int(Tab_connect[i,1]))+2]=k_global[2*(int(Tab_connect[i,0])):2*(int(Tab_connect[i,0]))+2,2*(int(Tab_connect[i,1])):2*(int(Tab_connect[i,1]))+2]+Mat_elem[i][0:2,2:4]
        k_global[2*(int(Tab_connect[i,1])):2*(int(Tab_connect[i,1]))+2,2*(int(Tab_connect[i,0])):2*(int(Tab_connect[i,0]))+2]=k_global[2*(int(Tab_connect[i,1])):2*(int(Tab_connect[i,1]))+2,2*(int(Tab_connect[i,0])):2*(int(Tab_connect[i,0]))+2]+Mat_elem[i][2:4,0:2]
    #k_reduit = np.delete(k_global, ddl_supprimer, axis=0)
    #k_reduit = np.delete(k_reduit, ddl_supprimer, axis=1)
    tous_ddl = np.arange(nombre_total_de_noeuds*2)
    ddl_present = np.setdiff1d(tous_ddl, ddl_supprimer)
    k_reduit = k_global[np.ix_(ddl_present, ddl_present)]
    return k_global, k_reduit

# Cette fonction assemble le vecteur force
def Assemblage_force(nombre_total_de_noeuds, Tab_connect, ddl_supprimer, force_elm):
    F_global = np.zeros((nombre_total_de_noeuds*2, 1))
    for i in range(np.shape(Tab_connect)[0]):
        F_global[2*(int(Tab_connect[i,0])):2*(int(Tab_connect[i,0]))+2, 0] += force_elm[i][0:2, 0]
        F_global[2*(int(Tab_connect[i,1])):2*(int(Tab_connect[i,1]))+2, 0] += force_elm[i][2:4, 0]
    #F_global = F_global.reshape((nombre_total_de_noeuds*2, 1))
    #F_reduit = np.delete(F_global, ddl_supprimer, axis=0)
    tous_ddl = np.arange(nombre_total_de_noeuds*2)
    ddl_present = np.setdiff1d(tous_ddl, ddl_supprimer)
    F_reduit = F_global[ddl_present, :]
    return F_global, F_reduit

# Cette fonction supprime les ddl du vecteur
def Spprime_ddl_vect(vecteur_colone, ddl_supprimer, ddl_present_non_bloque):
    #vecteur_reduit = np.delete(vecteur_colone, ddl_supprimer, axis=0)
    vecteur_reduit = vecteur_colone[ddl_present_non_bloque, :]
    return vecteur_reduit

# Cette fonction fait la repartition nodale, et creer le vecteur force de chaque element
# charge_repartie est un scalaire, position_force et charge_ponctuelle sont des vecteurs (ex: np.array([x1,x2,x3]))
def Force_element(charge_repartie, charge_ponctuelle, position_force, Tab_connect, nombre_total_de_noeuds):
    numero_element, x = Localise_element(position_force, Tab_connect)
    force_elm = [np.zeros([4, 1]) for i in range(np.shape(Tab_connect)[0])]
    for i in range(np.shape(charge_ponctuelle)[0]):
        if numero_element[i]>=0 : 
            Nx = N(x[i], Tab_connect[numero_element[i], 2])
            p = charge_ponctuelle[i]
            #F = np.array([[p*Nx[0,0]], [p*Nx[1,0]], [p*Nx[2,0]], [p*Nx[3,0]]])
            F = p*Nx
            force_elm[numero_element[i]] += F
    
    for elm in range(np.shape(Tab_connect)[0]):
        l = Tab_connect[elm, 2]
        q = np.array([[charge_repartie*l/2], [charge_repartie*l**2/12], [charge_repartie*l/2], [-charge_repartie*l**2/12]])
        force_elm[elm] += q
    return force_elm

# Cette fonction calcul les vecteurs force elementaire residuelle liée a la dynamique du vehicule + elastique Fb_t
def Force_dynamique_elastique_elment(Tab_connect, Mat_elem, numero_element, ub_t, ddl_global, fs_t, kbc):
    Fb_t_elm = []
    ue_t = np.zeros([4, 1])
    for e in range(len(Mat_elem)):
        if e == numero_element:
            Fb_e = fs_t + kbc @ ub_t    # Element sous la roue
        else:
            noeud_dep = int(Tab_connect[e, 0])
            noeud_arr = int(Tab_connect[e, 1])
            #ue_t[0:2] = ddl_pont[0][j][2*noeud_dep : 2*noeud_dep+2] 
            #ue_t[2:4] = ddl_pont[0][j][2*noeud_arr : 2*noeud_arr+2]
            ue_t[0:2] = ddl_global[2*noeud_dep : 2*noeud_dep+2] 
            ue_t[2:4] = ddl_global[2*noeud_arr : 2*noeud_arr+2]
            Fb_e = Mat_elem[e] @ ue_t
        Fb_t_elm.append(Fb_e)
    return Fb_t_elm

# Cette fonction permet de definir les charges ponctuelles et repartie qui sont sur la poutre
def Definition_charge():
    nbr_force = int(input(" Entrez le nombre de charge ponctuelle statique sur la poutre: "))
    input(" L'origine du repere est l'extremité gauche de la poutre. une charge est positive si elle est orienté vers le haut et negative dans le cas contraire ")
    force = np.zeros([nbr_force])
    position_force = np.zeros([nbr_force])
    for i in range(nbr_force):
        exposant = ("ère")*(i==0)+("ème")*(i!=0)
        force[i] = float(input(" Entrer la valeure de la "+str(i+1)+" "+str(exposant)+" force (en kN): F"+str(i+1)+" = "))*1000
        position_force[i] = float(input(" Entrer la position de la "+str(i+1)+" "+str(exposant)+" force (en m): x"+str(i+1)+" = "))
    q = float(input("Entrer la valeur de la charge uniformement repartie sur la poutre (en kN/m): q = "))*1000
    return force, position_force, q

# Cette fonction calcul les increments de deplacement du vehicule ainsi que la force qui agit sur le pont
# z = np.array([[z1],[z2]]), zp=vitesse et zpp=acceleration
def Deplacement_vehicule_force_agissant(z, zp, zpp, ub, rc, Nc):
    #global p, mw, Mv, kv, cv, kB, v, a0, a1, a2, a3, a4, a5, a6, a7
    z1, z1p, z1pp = z[0,0], zp[0,0], zpp[0,0]
    z2p, z2pp = zp[1,0], zpp[1,0]
    qs_t = np.array([[kv+kB, -kv],
                   [-kv,    kv]]) @ z
    qe_t = np.array([[-mw*(a2*z1p+a3*z1pp)-cv*(a4*(z1p-z2p)+a5*(z1pp-z2pp))],
                     [-Mv*(a2*z2p+a3*z2pp)-cv*(a4*(z2p-z1p)+a5*(z2pp-z1pp))]])
    qs1_t, qe1_t = qs_t[0,0], qe_t[0,0]
    qs2_t, qe2_t = qs_t[1,0], qe_t[1,0]
    qec_t = qe1_t+qe2_t
    qsc_t = qs1_t+qs2_t
    D = np.linalg.det(np.array([[kv+kB+a0*mw+a1*cv,        -kv-a1*cv],
                                [-kv-a1*cv,           kv+a0*Mv+a1*cv]]))
    delta_z = (-1/D)*np.array([[kv+a0*Mv+a1*cv],
                               [kv+a1*cv]]) * (p+kB*rc+kB*Nc.T @ ub)-(1/D)*np.array([[(qs1_t+qe1_t)*a0*Mv + (qsc_t+qec_t)*(kv+a1*cv)],
                                                                                           [(qs1_t+qe1_t)*(a0*mw+kB) + (qsc_t+qec_t)*(kv+a1*cv)]])
    kbc_vehicule = kB*(a0/D)*((Mv+mw)*(kv+a1*cv) + a0*Mv*mw)*Nc @ Nc.T  # kbc_vehicule = kbc - kb
    ps_t_plus_dt = -kB*(rc-((p+kB*rc)*(1/D)*(kv+a0*Mv+a1*cv)))*Nc 
    fs_t = kB*((1/D)*((qs1_t+qe1_t)*a0*Mv + (qsc_t+qec_t)*(kv+a1*cv))-z1)*Nc                                                                                      
    return delta_z, kbc_vehicule, ps_t_plus_dt, fs_t

# Cette fonction calcul la vitesse, le deplacement, l'acceleration pour une mise à jour à chaque pas de temps
def Vit_depl_accel(delta_u, u_t, up_t, upp_t):
    u_t_plus_dt = u_t + delta_u
    upp_t_plus_dt = a0*delta_u - a2*up_t - a3*upp_t   
    up_t_plus_dt = up_t + a6*upp_t + a7*upp_t_plus_dt  
    return u_t_plus_dt, up_t_plus_dt, upp_t_plus_dt,              

# Cette fonction calcul la matrice d'amortissement en supposant qu'elle est du type Rayleigh
def Amortissement(Mb, Kb):
    w_quare, phi = eigh(Kb, Mb) # Recuperation des pulsations et de la matrice modale
    if np.any(w_quare[:2]) <= 0 :
        raise np.linalg.LinAlgError("Les deux premieres valeurs propres doivent etre positive")
    if np.isclose(w_quare[0], w_quare[1]):
        raise ValueError(" Les deux frequences utilisés pour Rayleigh sont trop proche")
    w1 = np.sqrt(w_quare[0])
    w2 = np.sqrt(w_quare[1])                  
    alpha_0 = (2*w1*w2*(epsilon_2*w1-epsilon_1*w2))/(w1**2 - w2**2)
    alpha_1 = (2*(epsilon_1*w1-epsilon_2*w2))/(w1**2 - w2**2)
    Cb = alpha_0*Mb + alpha_1*Kb
    return Cb

# Cette fonction verifie si la matrice de rigidité est singuliere et arrete le programme car structure instable
def Matrice_siguliere(Kb_red):
    rang = np.linalg.matrix_rank(Kb_red)
    if rang < Kb_red.shape[0] :
        print("\n Matrice de rigité singuliere: Probleme de stabilité ou de modelisation. verifiez les appuis")
        raise np.linalg.LinAlgError("Le Programme va devoir s'arreter pour ne pas generer des resultats eronés.")
        # Dans mon futur programme, je dois separer les fontion analyse statique et dynamique, arreter l'analyse statique si Kb est singuliere mais poursuivre avec l'analyse dynamique(sauf si Kbc, matrice de rigité efective due à la presence du vehicule est singuliere)

# Cette fonction control les entées utilisateur
def Control_entree(v, nbr_elements, nbr_dt, tol, max_it):
    if v <= 0 :
        raise ValueError("La vitesse doit etre strictement positive")
    if nbr_elements < 1 : 
        raise ValueError("Il faut au moins un element par travée")
    if nbr_dt < 1 :
        raise ValueError(" Le nombre de Pas doit etre superieure ou egale à 1")
    if tol <= 0:
        raise ValueError("La tolerance doit etre positive")
    if max_it < 1 :
        raise ValueError("Le nombre maximal d'iteration doit etre superieure ou egale à 1")
        
# Cette fonction fait l'analyse dynamique du tablier modelisé par VBI
def Analyse_statique_dynamique():
    global p, mw, Mv, kv, cv, kB, v, a0, a1, a2, a3, a4, a5, a6, a7, epsilon_1, epsilon_2
    ddl_pont, ddl_vehicule = [[], [], []], [[], [], []] # Chaque liste contient les donnees relatives aux ddl du pont ou vehicule et les 03 sous listes de chaque liste contient les deplacement, vitesses et acceleration de tous les ddl
    U_glo_stat_plus_dynamique = []
    U_glo_stat_plus_vehicule = []       # contien les deplacement du pont sous ses charges (poids propre + ...) + celles du vehicule (statique) indispensable pour les lignes d'influences
    U_glo_vehicule_LI = []       # contien les deplacement du pont sous ses charges (poids propre + ...) + celles du vehicule (statique) indispensable pour les lignes d'influences
    position_vehicule = []  # Liste contenant l'historique des positions du vehicule
    force_de_contact = []
    historique_temporelle = []
    force_residulle_dynamique_plus_elastique = [] # Liste devrant contenir les forces ressiduelles liées à la dynamique de la suspension du vehicule + les forces elastiques internes de l'element de pont
    mw, Mv, kv, cv, kB, v, nbr_elements, nbre_dt, epsilon_1, epsilon_2, tol, max_it = Donnee_entree()
    Control_entree(v, nbr_elements, nbre_dt, tol, max_it)
    p = -g*(Mv + mw)
    Tab_connect, Mat_elem, Masse_elm, nombre_total_de_noeuds, ddl_intermediaire_supprimé = Definition_pout(nbr_elements)
    ddl_rive_supprimer = Numeros_DDL_a_supprimer(nombre_total_de_noeuds)
    ddl_supprimer = ddl_intermediaire_supprimé + ddl_rive_supprimer # Concatenation des deux listes
    Fp_stat, x_Fp_stat, q = Definition_charge()
    force_elm = Force_element(q, Fp_stat, x_Fp_stat, Tab_connect, nombre_total_de_noeuds)
    Mb_glo, Mb_red = Assemblage_matrice(nombre_total_de_noeuds, Tab_connect, Masse_elm, ddl_supprimer) # Matrice globale et reduite du pont
    Kb_glo_stat, Kb_red_stat = Assemblage_matrice(nombre_total_de_noeuds, Tab_connect, Mat_elem, ddl_supprimer)
    Matrice_siguliere(Kb_red_stat)
    
    Cb_red = Amortissement(Mb_red, Kb_red_stat)
    longueur_pont = Tab_connect[-1, 6]
    duree_analyse = longueur_pont/v  # Meme si la vitesse vari, je dois juste integrer l'acceleration du vehicule pour avoir la bonne formule qui calcule la durée d'analyse 
    subdivision_duree_analyse = np.linspace(0, duree_analyse, nbre_dt+1)
    dt = subdivision_duree_analyse[1] - subdivision_duree_analyse[0]
    a0 = 1/(beta*dt**2)
    a1 = gamma/(beta*dt)
    a2 = 1/(beta*dt)
    a3 = 1/(2*beta)-1
    a4 = gamma/beta - 1
    a5 = (dt/2)*((gamma/beta)-2)
    a6 = dt*(1-gamma)
    a7 = gamma*dt
    ub_t0, ubp_t0, ubpp_t0 = np.zeros([nombre_total_de_noeuds*2, 1]), np.zeros([nombre_total_de_noeuds*2, 1]), np.zeros([nombre_total_de_noeuds*2, 1]) # Initialisation des deplacement, vitesse et acceleration du pont (condition initiale)
    z_t0, zp_t0, zpp_t0 = np.zeros([2, 1]), np.zeros([2, 1]), np.zeros([2, 1]) # Initialisation des deplacement, vitesse et acceleration du vehicule (condition initiale)
    Fb_t = np.zeros([nombre_total_de_noeuds*2, 1])
    ddl_pont[0].append(ub_t0)
    ddl_pont[1].append(ubp_t0)
    ddl_pont[2].append(ubpp_t0)
    ddl_vehicule[0].append(z_t0)
    ddl_vehicule[1].append(zp_t0)
    ddl_vehicule[2].append(zpp_t0)
    force_residulle_dynamique_plus_elastique.append(Fb_t)
    position_vehicule.append(0)
    ub_t = np.zeros([4, 1])
    #u_stat = np.zeros([4, 1])
    #Fb_t_elm = [np.zeros([4,1]) for i in range(len(force_elm))]
    tous_ddl = np.arange(nombre_total_de_noeuds*2)
    ddl_present_non_bloque = np.setdiff1d(tous_ddl, ddl_supprimer)      # Retourne tous les numeros des noeuds DDL inconnu
    Ub_t_glo, Ubp_t_glo, Ubpp_t_glo = np.zeros([len(tous_ddl), 1]), np.zeros([len(tous_ddl), 1]), np.zeros([len(tous_ddl), 1])
    # Reponse du au force statiques sur le pont (poids propre et charge ponctuelle statique)
    Pb_t_plus_dt_stat_glo, Pb_t_plus_dt_stat_red = Assemblage_force(nombre_total_de_noeuds, Tab_connect, ddl_supprimer, force_elm)
    LU_Kb_red_stat = lu_factor(Kb_red_stat)        # factorisation de Kb_barre en LU
    U_red_stat = lu_solve(LU_Kb_red_stat, Pb_t_plus_dt_stat_red)   # inversion de la matrice de rigité du pont seul soumis au forces statique pour determiner les deplacement du pont
    U_glo_stat = np.zeros([len(tous_ddl), 1])
    U_glo_stat[ddl_present_non_bloque, :] = U_red_stat
    U_glo_stat_plus_dynamique.append(U_glo_stat.copy())
    U_glo_stat_plus_vehicule.append(U_glo_stat.copy())
    U_glo_vehicule_LI.append(np.zeros([len(tous_ddl), 1]))
    force_de_contact.append(-p)
    historique_temporelle.append(0)
    
    for j in range(len(subdivision_duree_analyse) - 1): # on parcourt le nombre de pas de temps
        xc = subdivision_duree_analyse[j+1]*v   # integrer l'acceleration deux fois pour avoir la formule adaptée au cas ou la vitesse n'est pas constante
        position_vehicule.append(xc)
        rc = r(xc)
        position_force = np.array([xc])
        numero_element, position_locale_de_la_roue = Localise_element(position_force, Tab_connect)
        if numero_element[0] >= 0 : # dans le cas ou on a plusieur roues, on enboite ce if dans une boucle qui parcour chaque roue
            Nc = N(position_locale_de_la_roue[0], Tab_connect[numero_element[0], 2])
            noeud_dep = int(Tab_connect[numero_element[0], 0])
            noeud_arr = int(Tab_connect[numero_element[0], 1])
            ub_t[0:2, :] = Ub_t_glo[2*noeud_dep : 2*noeud_dep+2, :] 
            ub_t[2:4, :] = Ub_t_glo[2*noeud_arr : 2*noeud_arr+2, :]
            #u_stat[0:2, :] = U_glo_stat[2*noeud_dep : 2*noeud_dep+2, :] 
            #u_stat[2:4, :] = U_glo_stat[2*noeud_arr : 2*noeud_arr+2, :]
            z = ddl_vehicule[0][j].copy()
            zp = ddl_vehicule[1][j].copy()
            zpp = ddl_vehicule[2][j].copy()
            delta_z, kbc_vehicule, ps_t_plus_dt, fs_t = Deplacement_vehicule_force_agissant(z, zp, zpp, ub_t, rc, Nc)
            Mat_elem_step = [K_elm.copy() for K_elm in Mat_elem]
            #kb = Mat_elem[numero_element[0]]
            kbc = Mat_elem_step[numero_element[0]] + kbc_vehicule
            Mat_elem_step[numero_element[0]] = kbc
            #pb_t_plus_dt = force_elm[numero_element[0]]
            #force_elm_step = [np.zeros([4, 1]) for i in range(np.shape(Tab_connect)[0])]
            #pb_t_plus_dt = force_elm_step[numero_element[0]]
            #force_elm_step[numero_element[0]] = pb_t_plus_dt + ps_t_plus_dt
            #force_elm_step[numero_element[0]] = ps_t_plus_dt
            force_elm_step = [F.copy() for F in force_elm]
            pb_t_plus_dt = force_elm_step[numero_element[0]]
            force_elm_step[numero_element[0]] = pb_t_plus_dt + ps_t_plus_dt

            Fb_t_elm = Force_dynamique_elastique_elment(Tab_connect, Mat_elem, numero_element[0], ub_t, Ub_t_glo, fs_t, kbc)
            #Fb_t_elm[numero_element[0]] = fs_t + kbc*ub_t
            Kb_glo, Kb_red = Assemblage_matrice(nombre_total_de_noeuds, Tab_connect, Mat_elem_step, ddl_supprimer)
            Pb_t_plus_dt_glo, Pb_t_plus_dt_red = Assemblage_force(nombre_total_de_noeuds, Tab_connect, ddl_supprimer, force_elm_step)
            Fb_t_glo, Fb_t_red = Assemblage_force(nombre_total_de_noeuds, Tab_connect, ddl_supprimer, Fb_t_elm)
            #force_residulle_dynamique_plus_elastique[j] = Fb_t_glo
            # les forces ci-dessous un un probleme d'indexation
            #delta_z, kbc_vehicule, ps_t_plus_dt, fs_t = Deplacement_vehicule_force_agissant(z_t0, zp_t0, zpp_t0, np.array([[0],[0],[0],[0]]), xc, Tab_connect[numero_element[0], 2])
            #Mat_elem[numero_element[0]] += kbc_vehicule 
            
            # Cette partie du code fait l'analyse statique du vehicule seul sur le pont 
            # charge_repartie est un scalaire, position_force et charge_ponctuelle sont des vecteurs (ex: np.array([x1,x2,x3]))
            force_elm_poid_vehicule = Force_element(0, np.array([p]), position_force, Tab_connect, nombre_total_de_noeuds)
            F_poid_vehicule_glo, F_poid_vehicule_red = Assemblage_force(nombre_total_de_noeuds, Tab_connect, ddl_supprimer, force_elm_poid_vehicule)
            U_red_stat_vehicule = lu_solve(LU_Kb_red_stat, F_poid_vehicule_red)   # inversion de la matrice de rigité du pont seul soumis au forces statique pour determiner les deplacement du pont
            U_glo_stat_vehicule = np.zeros([len(tous_ddl), 1])
            U_glo_stat_vehicule[ddl_present_non_bloque, :] = U_red_stat_vehicule
            U_glo_vehicule_LI.append(U_glo_stat_vehicule.copy())
            U_glo_stat_plus_vehicule.append(U_glo_stat.copy() + U_glo_stat_vehicule.copy())
        else:
            Nc = np.array([[0],[0],[0],[0]])
            ps_t_plus_dt = np.zeros([4, 1])
            fs_t = np.zeros([4, 1])
            Kb_glo, Kb_red = Assemblage_matrice(nombre_total_de_noeuds, Tab_connect, Mat_elem, ddl_supprimer)
            Pb_t_plus_dt_glo, Pb_t_plus_dt_red = Assemblage_force(nombre_total_de_noeuds, Tab_connect, ddl_supprimer, force_elm)
            Fb_t_glo = Kb_glo @ ddl_pont[0][j]
            Fb_t_red = Spprime_ddl_vect(Fb_t_glo, ddl_supprimer, ddl_present_non_bloque)
            U_glo_stat_plus_vehicule.append(U_glo_stat.copy())
            U_glo_vehicule_LI.append(np.zeros([len(tous_ddl), 1]))
        
        Kb_bar_t_plus_dt = a0*Mb_red + a1*Cb_red + Kb_red
        LU_Kb_bar_t_plus_dt = lu_factor(Kb_bar_t_plus_dt)        # factorisation de Kb_barre en LU
        Matrice_siguliere(Kb_bar_t_plus_dt)  # si la matrice est singuliere, le programme s'arrete
        
        # Fb_tl = force_residulle_dynamique_plus_elastique[j]
        Fb_t_red = Spprime_ddl_vect(force_residulle_dynamique_plus_elastique[j], ddl_supprimer, ddl_present_non_bloque)
        # Variable de l'instant t figer
        Ub_t_red_fige = Spprime_ddl_vect(ddl_pont[0][j].copy(), ddl_supprimer, ddl_present_non_bloque)
        Ubp_t_red_fige = Spprime_ddl_vect(ddl_pont[1][j].copy(), ddl_supprimer, ddl_present_non_bloque)
        Ubpp_t_red_fige = Spprime_ddl_vect(ddl_pont[2][j].copy(), ddl_supprimer, ddl_present_non_bloque)
        z_t_fige = ddl_vehicule[0][j].copy()
        zp_t_fige = ddl_vehicule[1][j].copy()
        zpp_t_fige = ddl_vehicule[2][j].copy() 
        #Ub_t_glo = ddl_pont[0][j].copy()
        
        # Variable iterative (celles qui vont converger)
        Ub_iter_red = Ub_t_red_fige.copy()
        Ubp_iter_red = Ubp_t_red_fige.copy()
        Ubpp_iter_red = Ubpp_t_red_fige.copy()
        z_iter = z_t_fige.copy()
        zp_iter = zp_t_fige.copy()
        zpp_iter = zpp_t_fige.copy()
        
        i = 1
        continuer = True
        while(continuer):
            Fb_bar_t_red = Fb_t_red - Mb_red @ (a2*Ubp_iter_red + a3*Ubpp_iter_red) - Cb_red @ (a4*Ubp_iter_red + a5*Ubpp_iter_red)  # Calcul avec les DDL reduit
            force_residuelle = Pb_t_plus_dt_red - Fb_bar_t_red
            norme_force_residuelle = np.linalg.norm(force_residuelle) # Norme de la force residuelle
            norme_force_reference = max(np.linalg.norm(Pb_t_plus_dt_red),1)
            norme_residuelle_relative = norme_force_residuelle/norme_force_reference
            
            # Extraction du deplacement local iteratif pour le contact
            Ub_t_glo[ddl_present_non_bloque, :] = Ub_iter_red
            if numero_element[0] >= 0:
                ub_t[0:2, :] = Ub_t_glo[2*noeud_dep : 2*noeud_dep+2, :] 
                ub_t[2:4, :] = Ub_t_glo[2*noeud_arr : 2*noeud_arr+2, :]
                fc = kB * ((Nc.T @ ub_t).item() + rc - z_iter[0,0])
            else:
                fc = 0
            
            if norme_residuelle_relative < tol and (fc >= 0 or numero_element[0] < 0):
                # --- CONVERGENCE ATTEINTE : On sauvegarde et on quitte le while ---
                Ubp_t_glo[ddl_present_non_bloque, :] = Ubp_iter_red
                Ubpp_t_glo[ddl_present_non_bloque, :] = Ubpp_iter_red
                ddl_pont[0].append(Ub_t_glo.copy())
                U_glo_stat_plus_dynamique.append(U_glo_stat.copy() + Ub_t_glo.copy())
                ddl_pont[1].append(Ubp_t_glo.copy())
                ddl_pont[2].append(Ubpp_t_glo.copy())
                ddl_vehicule[0].append(z_iter.copy())
                ddl_vehicule[1].append(zp_iter.copy())
                ddl_vehicule[2].append(zpp_iter.copy())
                force_residulle_dynamique_plus_elastique.append(Fb_t_glo.copy())
                force_de_contact.append(fc)
                historique_temporelle.append(subdivision_duree_analyse[j+1])
                continuer = False
            else:
                # --- CORRECTION ITÉRATIVE ---
                delta_U = lu_solve(LU_Kb_bar_t_plus_dt, force_residuelle)
                
                # Mise à jour Predictor-Corrector du Pont (Sans utiliser Vit_depl_accel !)
                Ub_iter_red += delta_U
                Ubp_iter_red += a1 * delta_U
                Ubpp_iter_red += a0 * delta_U
                
                if numero_element[0] >= 0:
                    # On ré-extrait le nouveau ub_t local après la correction
                    Ub_t_glo[ddl_present_non_bloque, :] = Ub_iter_red
                    ub_t[0:2, :] = Ub_t_glo[2*noeud_dep : 2*noeud_dep+2, :] 
                    ub_t[2:4, :] = Ub_t_glo[2*noeud_arr : 2*noeud_arr+2, :]
                    
                    # Éq 13 nécessite les variables FIGÉES au temps t
                    delta_z, kbc_vehicule, ps_t_plus_dt, fs_t = Deplacement_vehicule_force_agissant(z_t_fige, zp_t_fige, zpp_t_fige, ub_t, rc, Nc)
                    
                    # Mise à jour de l'état itératif du véhicule
                    z_iter = z_t_fige + delta_z
                    zpp_iter = a0*delta_z - a2*zp_t_fige - a3*zpp_t_fige
                    zp_iter = zp_t_fige + a6*zpp_t_fige + a7*zpp_iter
                    
                    # Recalcul des forces internes pour l'itération suivante
                    Fb_t_elm = Force_dynamique_elastique_elment(Tab_connect, Mat_elem, numero_element[0], ub_t, Ub_t_glo, fs_t, kbc)
                    Fb_t_glo, Fb_t_red = Assemblage_force(nombre_total_de_noeuds, Tab_connect, ddl_supprimer, Fb_t_elm)
                else:
                    Kb_glo, Kb_red = Assemblage_matrice(nombre_total_de_noeuds, Tab_connect, Mat_elem, ddl_supprimer)
                    Fb_t_glo = Kb_glo @ Ub_t_glo
                    Fb_t_red = Spprime_ddl_vect(Fb_t_glo, ddl_supprimer, ddl_present_non_bloque)
                    
                i += 1
                if i > max_it:
                    if fc < 0:
                        print("\n Perte de contact au pas "+str(j+1)+" à la position "+str(xc)+" lors de l'iteration "+str(i)+"")
                    else:
                        print("\n Non-convergence au pas de temps delta t = "+str(j+1)+" aprés "+str(max_it)+" iteration.")
                        print(f"Pas {j+1}, iteration {i}, résidu = {norme_force_residuelle: .3e}")
                    raise RuntimeError(f" Le programme va devoir s'arreter pour ne pas generer des resultats eronés")

    # ddl_pont contient la reponse dynamique du pont lors du passage du vehicule sans tenir compte des force statiques sur le pont            
    # ddl_vehicule contient les deplacements et acceleration verticaux de la caisse            
    # a supprimer U_glo_stat_plus_dynamique contient les deplacements totaux du à la dynamique du vehicule seul + les charges statique present sur le pont avant l'arrivé du vehicule
    # U_glo_stat_plus_vehicule contient les deplacements totaux du à la statique du vehicule seul (Ligne d'influence d'une charge ponctuelle se deplacant sur le pont) + les charges statique present sur le pont avant l'arrivé du vehicule
    # U_glo_vehicule_LI contient les deplacements totaux du à la statique du vehicule seul (Ligne d'influence d'une charge ponctuelle p = g*(Mv+mw) se deplacant sur le pont) 
    # position_vehicule Contient l'historique des positions occupées par le vehicule lors de sont passage sur le pont
    return ddl_pont, ddl_vehicule, U_glo_stat_plus_dynamique, U_glo_stat_plus_vehicule, U_glo_vehicule_LI, position_vehicule, force_de_contact, historique_temporelle


"""
            fc = kB*((Nc.T @ ub_t).item() + rc - z[0,0])*(numero_element[0] >= 0) + 0*(numero_element[0] < 0)
            if norme_residuelle_relative < tol and fc >= 0  :
                # sauvegarde des ddl
                Ubp_t_glo[ddl_present_non_bloque, :] = Ubp_t_red.copy()
                Ubpp_t_glo[ddl_present_non_bloque, :] = Ubpp_t_red.copy()
                ddl_pont[0].append(Ub_t_glo.copy())
                #U_glo_stat_plus_dynamique.append(U_glo_stat.copy() + Ub_t_glo.copy()) # Il contient les deplacement du pont sous charge statique + vehicule
                ddl_pont[1].append(Ubp_t_glo.copy())
                ddl_pont[2].append(Ubpp_t_glo.copy())
                ddl_vehicule[0].append(z.copy())
                ddl_vehicule[1].append(zp.copy())
                ddl_vehicule[2].append(zpp.copy())
                force_residulle_dynamique_plus_elastique.append(Fb_t_glo.copy())
                force_de_contact.append(fc)
                historique_temporelle.append(subdivision_duree_analyse[j+1])
                continuer = False
            else:
                delta_U = lu_solve(LU_Kb_bar_t_plus_dt, force_residuelle)   # inversion de la matrice de rigité equivalente pour determiner les increments de deplacement du pont
                Ub_t_glo[ddl_present_non_bloque, :], Ubp_t_glo[ddl_present_non_bloque, :], Ubpp_t_glo[ddl_present_non_bloque, :] = Vit_depl_accel(delta_U, Ub_t_red, Ubp_t_red, Ubpp_t_red)
                ddl_pont[0][j] = Ub_t_glo
                ddl_pont[1][j] = Ubp_t_glo
                ddl_pont[2][j] = Ubpp_t_glo
                Ub_t_red = Spprime_ddl_vect(Ub_t_glo, ddl_supprimer, ddl_present_non_bloque)
                Ubp_t_red = Spprime_ddl_vect(Ubp_t_glo, ddl_supprimer, ddl_present_non_bloque)
                Ubpp_t_red = Spprime_ddl_vect(Ubpp_t_glo, ddl_supprimer, ddl_present_non_bloque)
                if numero_element[0] >= 0 :
                    ub_t[0:2, :] = Ub_t_glo[2*noeud_dep : 2*noeud_dep+2, :] 
                    ub_t[2:4, :] = Ub_t_glo[2*noeud_arr : 2*noeud_arr+2, :]
                    delta_z, kbc_vehicule, ps_t_plus_dt, fs_t = Deplacement_vehicule_force_agissant(z_iter, zp_iter, zpp_iter, ub_t, rc, Nc)
                    z_new, zp_new, zpp_new = Vit_depl_accel(delta_z, z_t, zp_t, zpp_t)
                    z_iter = z_new
                    zp_iter = zp_new
                    zpp_iter = zpp_new
                    ddl_vehicule[0][j] = z
                    ddl_vehicule[1][j] = zp
                    ddl_vehicule[2][j] = zpp
                    Fb_t_elm = Force_dynamique_elastique_elment(Tab_connect, Mat_elem, numero_element[0], ub_t, Ub_t_glo, fs_t, kbc)
                    #Fb_t_elm[numero_element[0]] = fs_t + kbc*ub_t
                    Fb_t_glo, Fb_t_red = Assemblage_force(nombre_total_de_noeuds, Tab_connect, ddl_supprimer, Fb_t_elm)
                    #force_residulle_dynamique_plus_elastique[j] = Fb_t_glo
                else:
                    Kb_glo, Kb_red = Assemblage_matrice(nombre_total_de_noeuds, Tab_connect, Mat_elem, ddl_supprimer)
                    Fb_t_glo = Kb_glo @ Ub_t_glo
                    Fb_t_red = Spprime_ddl_vect(Fb_t_glo, ddl_supprimer, ddl_present_non_bloque)
                    #force_residulle_dynamique_plus_elastique[j] = Fb_t_glo
                i += 1
            if i > max_it :
                if fc < 0:
                    print("\n Perte de contact au pas "+str(j+1)+" à la position "+str(xc)+" lors de l'iteration "+str(i)+"")
                else:
                    print("\n Non-convergence au pas de temps delta t = "+str(j+1)+" aprés "+str(max_it)+" iteration.")
                    print(f"Pas {j+1}, iteration {i}, résidu = {norme_force_residuelle: .3e}")
                raise RuntimeError(f" Le programme va devoir s'arreter pour ne pas generer des resultats eronés") """
"""                # sauvegarde des ddl
                Ubp_t_glo[ddl_present_non_bloque, :] = Ubp_t_red.copy()
                Ubpp_t_glo[ddl_present_non_bloque, :] = Ubpp_t_red.copy()
                ddl_pont[0].append(Ub_t_glo.copy())
                #U_glo_stat_plus_dynamique.append(U_glo_stat.copy() + Ub_t_glo.copy()) # Il contient les deplacement du pont sous charge statique + vehicule
                ddl_pont[1].append(Ubp_t_glo.copy())
                ddl_pont[2].append(Ubpp_t_glo.copy())
                ddl_vehicule[0].append(z.copy())
                ddl_vehicule[1].append(zp.copy())
                ddl_vehicule[2].append(zpp.copy())
                force_residulle_dynamique_plus_elastique.append(Fb_t_glo.copy())
                continuer = False """
        
ddl_pont, ddl_vehicule, U_glo_stat_plus_dynamique, U_glo_stat_plus_vehicule, U_glo_vehicule_LI, position_vehicule, force_de_contact, historique_temporelle = Analyse_statique_dynamique()

# Cette fonction permet de determiner les sollicitation MTN
def MTN():
    print("MTN")

# Cette fonction determine le coefficient d'amplification dynamiuqe
def Resultat_dynamique(ddl_vehicule, ddl_pont, position_vehicule, U_glo_vehicule_LI):
    z2pp = []   # Acceleration de la caisse
    Ub_max = [] # Liste contenant les Deplacement maximal pour chacune des positions occupées par la force 
    numero_ddl_Ub_max = []  # en une position donné du vehicule, on a un deplacement ou rotation maximel quelque part sur le pont. cette liste contient le numero des DDL max pour chaque position de la force 
    """for i in range(len(ddl_vehicule[2])):
        z2pp.append(ddl_vehicule[2][i][1])
    max_z2pp = max(z2pp)
    indice = z2pp.index(max_z2pp)
    position = position_vehicule[indice]"""
    z2pp = np.array([ddl_vehicule[2][i][1, 0] for i in range(len(ddl_vehicule[2]))])
    i_max = np.argmax(np.abs(z2pp))
    max_z2pp = z2pp[i_max]
    position = position_vehicule[i_max]
    print(f"Pour une vitesse horizontal V = {v} m/s, L' acceleration maximale de la caisse du vehicule est {max_z2pp}. Elle se produit à la position x = {position}")
    # Evaluation du confort de roulement

    """for i in range(len(ddl_pont[0])):
        Ub_max.append(max(ddl_pont[0][i]))
        numero_ddl_Ub_max.append(ddl_pont[0][i].index(Ub_max[i]))"""
    for i in range(len(ddl_pont[0])):
        u = np.asarray(ddl_pont[0][i]).ravel()
        indice = np.argmax(np.abs(u))
        Ub_max_i = u[indice]
        Ub_max.append(Ub_max_i)
        numero_ddl_Ub_max.append(indice)
    indice_t = np.argmax(np.abs(Ub_max))
    ddl_Ub_max = Ub_max[indice_t]       # Reponse dynamique maximale
    numero_contenant_ddl_max = numero_ddl_Ub_max[indice_t]      # avec ce numero, je vais pouvoir identifier si il s'agit d'une translation(si nombre paire) ou d'une rotation(si nombre impair) 
    if numero_contenant_ddl_max % 2 == 0 :
        message = "translation"
    else:
        message = "rotation"
    position_ddl_Ub_max = position_vehicule[indice_t]
    # Repose statique à la meme section ou le max de la reponse dynamique se produit
    Ub_stat_max = U_glo_vehicule_LI[indice_t][numero_contenant_ddl_max, 0]
    print(f"\n La reponse (en deplacement) maximal due à la dynamique du vehicule est {ddl_Ub_max}, il s'agit d'une {message}. Elle se produit en x = {position_ddl_Ub_max}. La reponse lors de l'analyse statique du vehicule dans cette meme section est {Ub_stat_max}")
    if np.isclose(Ub_stat_max, 0.0):
        raise ZeroDivisionError("\n La reponse statique de reference est nulle.")
    DAF = abs(ddl_Ub_max)/abs(Ub_stat_max)      # Coefficient d'amplification dynamique
    #facteur_impact = (ddl_Ub_max - Ub_stat_max)/Ub_stat_max
    facteur_impact = DAF - 1
    print(f" Le Facteur d'impact est donc: {facteur_impact}")


Resultat_dynamique(ddl_vehicule, ddl_pont, position_vehicule, U_glo_vehicule_LI)
"""
 dans mon futur programme je dois appeler la demander à l'utilisateur si il souhaite avoir les resultats (LI, MTN, ...)
 tenant compte de la dynamique du vehicule. pour le calcul du facteur d'impact, je dois appeler la fonction Analyse_statique_dynamique
 avec comme parametre les charge statique du pont tous nulle et pour les autre resultats (vehicule + charge) les paraùetre
 de cette fonction(force statique-poid propre,...) seront non nul
 
 
 """