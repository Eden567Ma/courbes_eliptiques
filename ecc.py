# Sorbonne Université LU3IN024 2023-2024
# TME 5 : Cryptographie à base de courbes elliptiques
#
# Etudiant.e 1 : OUAKED Massilva 21212519

from math import sqrt
import matplotlib.pyplot as plt
from random import randint

# Fonctions utiles

def exp(a, N, p):
    """Renvoie a**N % p par exponentiation rapide."""
    def binaire(N):
        L = list()
        while (N > 0):
            L.append(N % 2)
            N = N // 2
        L.reverse()
        return L
    res = 1
    for Ni in binaire(N):
        res = (res * res) % p
        if (Ni == 1):
            res = (res * a) % p
    return res


def factor(n):
    """ Return the list of couples (p, a_p) where p is a prime divisor of n and
    a_p is the p-adic valuation of n. """
    def factor_gen(n):
        j = 2
        while n > 1:
            for i in range(j, int(sqrt(n)) + 1):
                if n % i == 0:
                    n //= i
                    j = i
                    yield i
                    break
            else:
                if n > 1:
                    yield n
                    break

    factors_with_multiplicity = list(factor_gen(n))
    factors_set = set(factors_with_multiplicity)

    return [(p, factors_with_multiplicity.count(p)) for p in factors_set]


def inv_mod(x, p):
    """Renvoie l'inverse de x modulo p."""
    return exp(x, p-2, p)


def racine_carree(a, p):
    """Renvoie une racine carrée de a mod p si p = 3 mod 4."""
    assert p % 4 == 3, "erreur: p != 3 mod 4"

    return exp(a, (p + 1) // 4, p)

def est_premier(p) :
    for i in range(2, p-1) :
        if ( p % i ) == 0 :
            return False
    return True

    
# Fonctions demandées dans le TME

def est_elliptique(E):
    """
    Renvoie True si la courbe E est elliptique et False sinon.

    E : un triplet (p, a, b) représentant la courbe d'équation
    y^2 = x^3 + ax + b sur F_p, p > 3
    """
    (p, a,b)=E 
    delta=(4*a**3 +27*b**2)%p 
    if delta ==0 :
        return False
    n= p+1
    for x in range (p):
        w =(x**3 + a*x +b)% p 
        
    for y in range (2,int(n%0.5)+1):
        if n%y == 0 and p%y ==0:
            return False 

    return True

def pgcd(x,y) :
    while a!=b: 
        d=abs(b-a) 
        b=a 
        a=d 
    return d

def point_sur_courbe(P, E):
    """Renvoie True si le point P appartient à la courbe E et False sinon.""" 
    
    
    if P is None:
        # None represents the point at infinity.
        return True
    print("Courbes")
    x, y = P
    p, a, b = E
    
    return (y * y - x * x * x - a * x - b) % p == 0



def symbole_legendre(a, p):
    """Renvoie le symbole de Legendre de a mod p."""
    return pow(a, (p - 1) // 2, p)


    
def cardinal(E):
    """Renvoie le cardinal du groupe de points de la courbe E."""
    p, a, b = E
    cpt = 1 # l'infini
    for x in range(p):
            z = (x**3 + x * a + b) % p  #eviter la foction pow 

            legendre = symbole_legendre(z, p)
            
            if legendre == 1: cpt += 2           
            
            elif legendre == 0: cpt += 1
    #print(cpt)
    return cpt



    




def liste_points(E):
    """Renvoie la liste des points de la courbe E."""
    p, a, b = E
    resultat = []
    y1 = []

    for y in range((p//2)+1) :
        y1.append(pow(y,2,p))

    y1 = set(y1)

    for x in range(p):
        z = (pow(x, 3, p)+(a*x)%p+b)%p
        if z in y1 :
            if z != 0 : 
                y, others = pow(z, (p+1)//4, p), -pow(z, (p+1)//4, p)
                resultat.append((x, y))
                resultat.append((x, others%p))
            
            else :
                resultat.append((x,0))
        
    resultat.append(())
    return resultat


def cardinaux_courbes(p):
    """Calcule les cardinaux possibles pour une courbe elliptique définie sur Fp."""
    D = {}

    for b in range(p):
       
        for a in range(p):
            E = p, a, b
            
            if est_elliptique(E): 
               
                card = cardinal(E)
                if card in D :
                    D[card] += 1
                else:
                    D[card] = 1
    return D





def dessine_graphe(p):
    """Dessine le graphe de répartition des cardinaux des courbes elliptiques définies sur F_p."""
    bound = int(2 * sqrt(p))
    C = [c for c in range(p + 1 - bound, p + 1 + bound + 1)]
    D = cardinaux_courbes(p)

    plt.bar(C, [D[c] for c in C], color='b')
    plt.show()


def moins(P, p):
    """Retourne l'opposé du point P mod p."""
    x,y = P
    x = x % p
    y = y% p
    P2 = (-x, y)
    return P2


def est_egal(P1, P2, p):
    """Teste l'égalité de deux points mod p."""
    if (P1 == () and P2 == ()) :
        return True
    elif ( P1 == () and P2 != () ) :
        return False 
    elif (P1 != () and P2 ==()) :
        return False 
    else :
        x1, y1 = P1
        x2, y2 = P2
        first_egalite = x1 % p
        second_egalite = x2 % p
        fisrt1_egalite = y1 % p
        secondd_egalite = y2 % p
        return ( first_egalite == second_egalite) and (fisrt1_egalite == secondd_egalite)


def est_zero(P):
    """Teste si un point est égal au point à l'infini."""
    return (P==())


def addition(P1, P2, E):
    """Renvoie P1 + P2 sur la courbe E."""
    '''p, a, b = E
    if est_zero(P1):
        return P2
    if est_zero(P2) :
        return P1
    #si P1 est l'opposé de P2, P1+P2 = O
    if est_egal(moins(P1, p),P2,p) :
        return ()
    x1, y1 = P1
    x2, y2 = P2
    if est_egal(P1,P2,p) :
        la = (3 *pow(x1,2,p) + a)%p
        la = la*inv_mod((2*y1)%p,p)
    #sinon les deux points sont differents
    else :
        la = (y2-y1)%p
        la = la * inv_mod((x2-x1)%p,p)
    x3 = (pow(la,2,p) - x1 -x2) % p
    y3 = (la*(x1-x3)-y1) % p
    return x3,y3'''





def multiplication_scalaire(k, P, E):
    """Renvoie la multiplication scalaire k*P sur la courbe E."""

    p, a, b = E
    c = bin(k)
    resultat = P
    if k < 0 :
        uP = moins(P,p)
        k = -k
    else :
        up = P

    kP = ()
    P2 = P 
    L = L = bin(k >> 1)
    for i in range(len(L)) :
        if L[-1-i] == b :
            break
        P2 = addition(P2,P2,E)
        if L[-1-i] == 1 :
            kP = addition(kP, P2, E)
    if k % 2 == 1 : 
        kP = addition(kP, P, E)
    return kP