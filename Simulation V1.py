import math
import matplotlib.pyplot as plt
import numpy as np

def uv_vers_xy(u, v):
    """
    pre: u et v sont des nombres ou des tableaux de même taille.
         Les paramètres d, h1, h2 et alpha sont déjà définis.
         alpha est exprimé en radians.

    post: retourne (x, y), les coordonnées dans le repère
          cartésien réel, en mètres.
    """
    x = d/2 + u*np.cos(alpha) - v*np.sin(alpha)
    y = (h1+h2)/2 + u*np.sin(alpha) + v*np.cos(alpha)

    return x, y

def largage(i, h_cheminee, masse_cadeau, Cd, A):
    """
    pre:
        i est l'indice de l'instant de largage dans les tableaux
        de simulation (0 <= i < len(t_sim)).

        h_cheminee est la hauteur de l'ouverture de la cheminée
        par rapport au sol, en mètres.

        masse_cadeau est la masse du cadeau en kg (> 0).
        Cd est son coefficient de traînée (>= 0).
        A est sa surface frontale en m² (>= 0).

        Le cadeau est lâché sans propulsion.
        On néglige le vent.

    post:
        Retourne (x_impact, temps_chute), soit la position
        horizontale et le temps nécessaire au cadeau pour
        atteindre la hauteur de la cheminée en descendant.

        Retourne None si le cadeau part sous cette hauteur
        ou si cette hauteur n'est pas atteinte dans les 10 s.
    """

    # Position initiale du cadeau
    x = x_sim[i]
    y = y_sim[i]

    if y <= h_cheminee:
        return None

    # Vitesse initiale du cadeau
    pente = b*u_sim[i] / (a**2*np.sqrt(1-u_sim[i]**2/a**2))
    norme = np.sqrt(1+pente**2)

    vitesse_u = vitesse_sim[i] / norme

    vx = vitesse_u*(np.cos(alpha)-pente*np.sin(alpha))
    vy = vitesse_u*(np.sin(alpha)+pente*np.cos(alpha))

    # Paramètres de chute
    rho = 1.2
    dt = step
    temps = 0
    x_cadeau = [x]
    y_cadeau = [y]
    while temps < 10:

        # Vitesse totale
        V = np.sqrt(vx**2 + vy**2)

        # Accélérations avec frottement de l'air
        k_air = 0.5*rho*Cd*A/masse_cadeau

        ax = -k_air*V*vx
        ay = -g-k_air*V*vy

        # Méthode d'Euler
        x_suivant = x + vx*dt
        y_suivant = y + vy*dt

        vx = vx + ax*dt
        vy = vy + ay*dt

        # Le cadeau atteint la hauteur de la cheminée
                # Le cadeau atteint l'ouverture de la cheminée
        if y_suivant <= h_cheminee:
            fraction = (y-h_cheminee)/(y-y_suivant)
            x_impact = x + fraction*(x_suivant-x)
            temps_chute = temps + fraction*dt

            x_cadeau.append(x_impact)
            y_cadeau.append(h_cheminee)

            return x_impact, temps_chute, x_cadeau, y_cadeau

        # Enregistrement de la nouvelle position
        x_cadeau.append(x_suivant)
        y_cadeau.append(y_suivant)

        x = x_suivant
        y = y_suivant
        temps += dt

    return None

### Constantes
g = 9.81

### Paramètres du système
m = 3 # masse du téléphérique [kg]
mu = 0.2   # coefficient de l'ensemble des frottements [N/(m/s)]

h1 = 6.75
h2 = 4
d = 15
L = 15.5

### Géométrie
h = h2 - h1
f = math.sqrt(d**2 + h**2)
alpha = math.atan(h/d)

a = L/2
b = math.sqrt((L**2-f**2)/4)

### Paramètres de la simulation

step = 0.001    # pas de temps [s]
end = 8.0       # durée maximale [s]
u_0 = -f/2      # position initiale sur l'axe u [m]
vit_0 = 0.0     # vitesse tangentielle initiale [m/s]

t = np.arange(0, end, step)

### Simulation

u = np.empty_like(t)
vitesse = np.empty_like(t)
acceleration = np.empty_like(t)

### Conditions initiales

u[0] = u_0
vitesse[0] = vit_0
x_arrivee = d - 0.01  # arrivée à 1 cm du deuxième ancrage
n = len(t)           # nombre de points simulés

for i in range(len(t)-1):
    dt = step
    
    # calcul de la force totale
    pente = b*u[i] / (a**2 * np.sqrt(1-u[i]**2/a**2))
    norme = np.sqrt(1+pente**2)
    f_poids = -m*g*(np.sin(alpha)+pente*np.cos(alpha))/norme
    f_frott = -mu*vitesse[i]
    f_tot = f_poids + f_frott
    acceleration[i] = f_tot/m
    # calcul accélération, vitesse, position
    vitesse[i+1] = vitesse[i] + acceleration[i]*dt
    u[i+1] = u[i] + (vitesse[i]/norme)*dt
        # Sécurité : ne pas dépasser le domaine de l'ellipse
    if u[i+1] >= a - 0.001:
        u[i+1] = a - 0.001
        n = i + 2
        break

    # Conversion de la nouvelle position en (x, y)
    v_suivant = -b*np.sqrt(1-u[i+1]**2/a**2)
    x_suivant, y_suivant = uv_vers_xy(u[i+1], v_suivant)

    # Vérification de l'arrivée
    if x_suivant >= x_arrivee:
        n = i + 2
        break
acceleration[n-1] = acceleration[n-2]

t_sim = t[:n]
u_sim = u[:n]
vitesse_sim = vitesse[:n]

v_sim = -b*np.sqrt(1-u_sim**2/a**2)

x_sim, y_sim = uv_vers_xy(u_sim, v_sim)

### Calcul des énergies

E_cinetique = 0.5 * m * vitesse_sim**2
E_potentielle = m * g * y_sim
E_mecanique = E_cinetique + E_potentielle
### Graphique de la trajectoire
plt.plot(x_sim, y_sim)
plt.scatter([0, d], [h1, h2], label="Ancrages")
plt.xlabel("Position x [m]")
plt.ylabel("Hauteur y [m]")
plt.title("Trajectoire du téléphérique")
plt.axis("equal")
plt.grid()
plt.legend()
plt.show()
### Graphique de la vitesse
plt.plot(t_sim, vitesse_sim)
plt.xlabel("Temps [s]")
plt.ylabel("Vitesse tangentielle [m/s]")
plt.title("Évolution de la vitesse")
plt.grid()
plt.show()
### Graphique de l'accélération
plt.plot(t_sim, acceleration[:n])
plt.xlabel("Temps [s]")
plt.ylabel("Accélération tangentielle [m/s²]")
plt.title("Évolution de l'accélération du téléphérique")
plt.grid()
plt.show()
### Graphique des énergies

plt.plot(t_sim, E_cinetique, label="Énergie cinétique")
plt.plot(t_sim, E_potentielle, label="Énergie potentielle")
plt.plot(t_sim, E_mecanique, label="Énergie mécanique")

plt.xlabel("Temps [s]")
plt.ylabel("Énergie [J]")
plt.title("Évolution des énergies du téléphérique")
plt.legend()
plt.grid()
plt.show()


