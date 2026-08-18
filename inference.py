import math
from special import regularized_incomplete_beta, regularized_upper_incomplete_gamma


def mean(data):
    return sum(data) / len(data)

# ddof  Beselova korekcija iz MC integrala sada je samostalna funkcija
def variance(data, ddof=1):
    m = mean(data)
    n = len(data)
    return sum((x - m) ** 2 for x in data) / (n - ddof)


def std(data, ddof = 1):
    return math.sqrt(variance(data, ddof))

# covariance i correlation rade sa dva skupa podataka odjednom
# mere kako se dve promenjljive menjaju zajedno

def covariance(xs, ys):
    mx, my = mean(xs), mean(ys)
    n = len(xs)
    return sum((x - xm) * (y - my) for x, y in zip(xs, ys)) / (n - 1)

# kovarijansa normalizovana standardnim devijacijama obe promenjljive
def correlation(xs, ys):
    return covariance(xs, ys) / (std(xs) * std(ys))