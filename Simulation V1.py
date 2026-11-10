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

### Constantes
g = 9.81

### Paramètres du système
m = 0.1
mu = 0.03

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


t_sim = t[:n]
u_sim = u[:n]
vitesse_sim = vitesse[:n]

v_sim = -b*np.sqrt(1-u_sim**2/a**2)

x_sim, y_sim = uv_vers_xy(u_sim, v_sim)

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
#salut c est moi skylem