"""
    special.py
    ------------------------
    Regularizovana nepotpuna beta funkcija i nepotpuna gama funkcija

"""

import math
#
# nepotpuna beta funkcija
#
# neprekidni razlomak
def _betacf(a,b,x):
    MAXIT = 200    # maksimalan broj spratova neprekidnog razlomka
    EPS = 3e-9     # epsilon- prag preciznosti
    FPMIN = 1e-30   # zastida od deljenja sa nulom

    qab = a+b
    qap = a + 1.0
    qam = a - 1.0

    # pocetne vrednosti za Lencov alogoritam
    # dve pomocne promenjljive koje se iterativno azuriraju
    # svaki -sprat- neprekidnog razlomka menja i c i d
    # obe su zajeno trenutni razlomak izracunat do tog -sprata-
    c = 1.0
    d = 1.0 - qab * x/qap
    if abs(d) < FPMIN:   # zastita deljenja sa nulom
        d = FPMIN

    d = 1.0 / d   # Lorenc
    h = d     # akomulira konacan rezultat kroz petlju -
              # pocinje sa istom vrednoscu kao d pa se mnozi kroz svaki sprat

    # Lorencov alogoritam za beta funkciju
    for m in range(1, MAXIT +1):
        m2 = 2 * m

        # prvi korak razlomka paran
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c

        # drugi korak razlomka neparni
        aa = -(a +m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta

        if abs(delta - 1.0) < EPS:   #ako nova promena nista ne menja- prekid
            break

    return h

def regularized_incomplete_beta (a, b ,x):

   # granicni slucajevi
   if x <= 0.0:
        return 0.0
   if x >= 1.0:
        return 1.0
    # predfaktor racunat preko alogoritma
   ln_beta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
               + a * math.log(x)
               + b * math.log(1.0 - x))
   front = math.exp(ln_beta)

   # bira se grana koja brze konvergira
   # za pogresnu stran u- identittet ->  I_x(a,b) = 1 - I_(1-x)(b,a)
   if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a,b,x) / a
   else:
        return 1.0 - front * _betacf(b, a, 1.0 - x) / b


# kvantilna funkcija
# inv. CDF Beta(a, b) raspodela preko binarne pretrage
# za verovatnocu p vraca x tako da CDF(x) = p
def beta_quantile(p, a, b, tol=1e-8):
    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if regularized_incomplete_beta(a, b, mid) < p:
            lo = mid      # pravi odgovor je desno od mid
        else:
            hi = mid      # pravi odgovor je levo ili tacno na mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2.0

#
# nepotpuna gama funkcija
#

def _gamma_series(a, x):
    # P(a,x) preko razvoja u red - konvergira brzo kada je x < a+1.
    ITMAX = 300    # broj koraka
    EPS = 3e-9   #prag preciznosti za rani prekid

    if x <= 0.0:
        return 0.0

    ap = a           # trenutni eksponent u redu
    total = 1.0 / a     # zbir svih clanova reda
    delta = total          # trenutni clan koji se dodaje
    for _ in range(ITMAX):
        ap +=1.0
        delta *= x / ap    # sledeci clan se racuna iz predhodnog
        total += delta
        if abs(delta) < abs(total) * EPS:  # rani prekid
            break

    # predfaktor i zavrsetak
    log_prefactor = -x + a * math.log(x) - math.lgamma(a)
    return total * math.exp(log_prefactor)

# neprekidni razlomak za gama funkciju
# Q(a, x) preko neprekidnog razlomka - konvergira brzo kada je x >= a+1
# Lencov princip -- druga formula
def _gamma_cf(a, x):
    ITMAX = 300
    EPS = 3e-9
    FPMIN = 1e-30

    b = x + 1.0 - a
    c = 1.0 / FPMIN
    d = 1.0 / b
    h = d

    for i in range(1, ITMAX + 1):
        an = -i * (i-a)
        b += 2.0
        d = an * d + b

        if abs(d) < FPMIN:
            d = FPMIN
        c = b + an / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            break

    log_prefactor = -x + a * math.log(x) - math.lgamma(a)
    return math.exp(log_prefactor) * h

# P(a, x) - donja nepotpuna gama funkcija ("udeo od 0 do x")
def regularized_lower_incomplete_gamma(a, x):
    if x < 0.0 or a <= 0.0:
        raise ValueError("a mora biti > 0, x mora biti >= 0")
    if x == 0.0:
        return 0.0

    if x < a + 1.0:
        return _gamma_series(a, x)        # red konvergira brze ovde
    else:
        return 1.0 - _gamma_cf(a, x)    #inace koristi razlomak -> pa oduzmi od 1

# Q(a, x) = 1 - P(a, x) - gornja nepotpuna gama funkcija.
# Koristi se DIREKTNO kao p-vrednost hi-kvadrat testa.
def regularized_upper_incomplete_gamma(a, x):
    if x < 0.0 or a <= 0.0:
        raise ValueError("a mora biti > 0, x mora biti >= 0")
    if x == 0.0:
        return 1.0

    if x < a + 1.0:
        return 1.0 - _gamma_series(a, x)    #oduzmi od 1 da dobijes Q
    else:
        return _gamma_cf(a, x)          # razlomak racuna Q direktno


