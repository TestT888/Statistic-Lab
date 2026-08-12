

import math
import random
from distributions import normal


# Monte Carlo procena broja pi
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

# Metropolis-Hastings
# Pravilo prihvatanja
# prihvati x' ako je: log(slucajan_broj_0_do_1)<log_density(x')-log_density(x)

def metropolis_hastings(log_density, initial, n_samples, proposal_std=1.0, burn_in=500):
    samples = []
    current = initial
    current_log_p = log_density(current)

    total_steps = n_samples + burn_in   #burn-in period

    for i in range(total_steps):            #provera da li je prosao burn-in period
        proposal = normal(current, proposal_std)
        proposal_log_p = log_density(proposal)

        log_accept_ratio = proposal_log_p - current_log_p  # pozitivno predlog je gusci, negativno predlog je redji
        # ako je predlog bolji ili jednako dobar  uvek se prihvata
        # ako je predlog losiji prihvati ga samo sa verovatnocom -> exp(log_acept_ratio)
        # sto je predlog losiji exp(ratio) je blize nuli
        if log_accept_ratio >= 0 or random.random() < math.exp(log_accept_ratio):
            current = proposal
            current_log_p = proposal_log_p
        if i>= burn_in:
            samples.append(current)        # tek posle burn-in perioda pozicije se cuvaju


    return samples


