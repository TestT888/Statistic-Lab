# Monte Carlo procena broja pi

import math
import random
from getopt import error


def monte_carlo_pi(n):
    inside = 0
    for _ in range(n):
        x = random.uniform(-1, 1)
        y = random.uniform(-1, 1)
        if x * x +y * y < 1.0:
            inside += 1

    p_hat = inside / n
    pi_estimate = 4* p_hat
    se = 4 * math.sqrt(p_hat * (1 - p_hat) /n)
    return pi_estimate, se


def monte_carlo_integral (f, a, b, n):
    samples = [f(random.uniform(a, b)) for _ in range(n)]

    # samples = []
    # for _ in range(n)
    #      x = random.uniform (a, b)
    #      samples.append(f(x))

    mean_f = sum(samples) / n

    # varijansa izracunatih vrednosti f(x)
    # koliko su vrednosti razbacane oko svog proseka
    # Beselova korekcija
    var_f = sum ((s-mean_f)**2 for s in samples)/(n-1)
    integral_estimate = (b-a)*mean_f
    se = (b-a)*math.sqrt(var_f/n)
    return integral_estimate, se

def clt_demo(sampler, sample_size, num_trials):
    means = []
    for _ in range(num_trials):
        sample = [sampler() for _ in range (sample_size)]
        means.append(sum(sample) / sample_size)
    return means

