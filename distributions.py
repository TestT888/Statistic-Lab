"""
distributions.py
-----------------
Generatori slucajnih uzoraka, izgradjeni iskljucivo preko `random`
(koji daje uniformnu raspodelu) i `math` modula.

Ideja: uniformna raspodela je "sirovina" - preko matematickih
transformacija (Box-Muller, inverzna transformacija...) iz nje
se prave sve ostale raspodele koje su potrebne.
"""
import random
import math

# Cuva "rezervni" broj iz Box-Muller transformacije izmedju poziva
# normal() funkcije - formula racuna DVA broja odjednom, ovde
# pamti drugi za sledeci poziv
_spare_gaussian = None

# Standardna normalna (Gausova) raspodela preko Box-Muller transformacije.
# mu = prosek, sigma = standardna devijacija
def normal(mu=0.0, sigma=1.0):
    global _spare_gaussian
    if _spare_gaussian is not None:
        z =  _spare_gaussian
        _spare_gaussian = None
        return mu+sigma*z

    u1 =random.random()
    u2 =random.random()
    u1 = max(u1, 1e-12)

    r = math.sqrt(-2.0 * math.log(u1))
    z0 = r * math.cos(2 * math.pi * u2)
    z1 = r * math.sin(2 * math.pi * u2)

    _spare_gaussian = z1
    return mu+sigma*z0


#generise uniformnu raspodelu
def uniform(a= 0.0, b=1.0):
    return a+(b-a)*random.random()

# inverzna transformacija
# eksponencijalna raspodela za modelovanje vremena cekanja
# veci- rate - vreme cekanja manje
def exponential(rate=1.0):
    u = random.random()
    u = max(u, 1e-12)
    return -math.log(u) / rate

# diskretne rspodele

# simulacija bacanja novcica
def bernoulli(p=0.5):
    return 1 if random.random() < p else 0
#simulacija n bacanja
def binomial(n, p=0.5):
    return sum(bernoulli(p) for _ in range(n))



