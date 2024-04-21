import numpy as np
from math import e
from scipy.stats import norm


d1 = lambda S, K, r, q, sigma, T: (np.log(S / K) + (r - q + (sigma**(2) * 0.5)) * T)/ \
                                    (sigma * np.sqrt(T))

d2 = lambda S, K, r, q, sigma, T: d1(S, K, r, q, sigma, T) - (sigma * np.sqrt(T))

price_call = lambda S, K, r, q, sigma, T: (S * e**(-q * T) * norm.cdf(d1(S, K, r, q, sigma, T)))\
                                        - (K * e**(-r * T) * norm.cdf(d2(S, K, r, q, sigma, T)))

price_put = lambda S, K, r, q, sigma, T: (K * e**(-r * T) * norm.cdf(-d2(S, K, r, q, sigma, T)))\
                                        - (S * e**(-q * T) * norm.cdf(-d1(S, K, r, q, sigma, T)))

delta_call = lambda S, K, r, q, sigma, T: e**(-q * T) * norm.cdf(d1(S, K, r, q, sigma, T))

delta_put = lambda S, K, r, q, sigma, T: e**(-q * T) * (norm.cdf(d1(S, K, r, q, sigma, T)) - 1)

gamma = lambda S, K, r, q, sigma, T: (e**(-q * T) / (S * sigma * np.sqrt(T))) * norm.pdf(d1(S, K, r, q, sigma, T))

vega = lambda S, K, r, q, sigma, T: S * e**(-q * T) * np.sqrt(T) * norm.pdf(d1(S, K, r, q, sigma, T))

theta_call = lambda S, K, r, q, sigma, T: ((-(S * sigma * e**(-q * T)) / (2 * np.sqrt(T))) * norm.pdf(d1(S, K, r, q, sigma, T)))\
                                        - (r * K * e**(-r * T) * norm.cdf(d2(S, K, r, q, sigma, T)))\
                                        + (q * S * e**(-q * T) * norm.cdf(d1(S, K, r, q, sigma, T)))

theta_put = lambda S, K, r, q, sigma, T: ((-(S * sigma * e**(-q * T)) / (2 * np.sqrt(T))) * norm.pdf(d1(S, K, r, q, sigma, T)))\
                                        + r * K * e**(-r * T) * norm.cdf(-d2(S, K, r, q, sigma, T))\
                                        - q * S * e**(-q * T) * norm.cdf(-d1(S, K, r, q, sigma, T))

rho_call = lambda S, K, r, q, sigma, T: K * T * e**(-r * T) * norm.cdf(d2(S, K, r, q, sigma, T))

rho_put = lambda S, K, r, q, sigma, T: -K * T * e**(-r * T) * norm.cdf(-d2(S, K, r, q, sigma, T))