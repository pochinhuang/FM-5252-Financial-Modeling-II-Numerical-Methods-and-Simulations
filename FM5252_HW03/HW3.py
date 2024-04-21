import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize

class bs:

    """
    Part One:
        Black Scholes, bisection, and Newtons method.
    
    Part Two:
        Define a error function to find the mean of squared distances between ws and implied volitilities.
        Pass this error function to scipy minimizer, it minimized the error with a, b, ρ, m, σ.
        Plot it with a, b, ρ, m, σ.

    
    S: Stock price
    K: Strike price
    r: Risk free rate
    q: Dividend
    sigma: Volitility
    T: Time until option expiration
    option_type: Call or Put

    
    References:
    https://www.cnblogs.com/sljsz/p/15664358.html
    https://0809zheng.github.io/2021/08/23/minimize.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
    https://stackoverflow.com/questions/55187816/scipy-optimize-constraints

    """

    def __init__(self, S, K, r, q, sigma, T, option_type):
        self.S = S
        self.K = K
        self.r = r
        self.q = q
        self.sigma = sigma
        self.T = T
        self.option_type = option_type
        self.d1 = (np.log(S / K) + (r - q + (sigma**(2) * 0.5)) * T)/ (sigma * np.sqrt(T))
        self.d2 = self.d1 - (sigma * np.sqrt(T))
 

 # Black–Scholes model
    def price(self):
        if self.option_type == "call":
            price = (self.S * np.exp(-self.q * self.T) * norm.cdf(self.d1))\
                    - (self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2))
            return price
            
        elif self.option_type == "put":
            price = (self.K * np.exp(-self.r * self.T) * norm.cdf(-self.d2))\
                    - (self.S * np.exp(-self.q * self.T) * norm.cdf(-self.d1))
            return price
        
        else:
            raise ValueError("Invalid option type. Use 'call' or 'put'.")   


    # calculate vega
    def vega(self):
        vega = self.S * np.exp(-self.q * self.T) * np.sqrt(self.T) * norm.pdf(self.d1)
        return vega


    # bisection method
    def bisection(self, option_price, high = 1.0, low = 0.001, tolerance = 0.00000001):
  
        while high - low > tolerance:

            mid = (low + high) / 2

            option = bs(self.S, self.K, self.r, self.q, mid, self.T, self.option_type)
            mid_price = option.price()
            diff = option_price - mid_price

            if diff > 0:
                low = mid
            else:
                high = mid

        return mid
    

    # Newton’s method
    def newton(self, option_price, guess = 0.5, diff = 1, tolerance = 0.00000001):
        
        while abs(diff) >= tolerance:

            option = bs(self.S, self.K, self.r, self.q, guess, self.T, self.option_type)
            price = option.price()
            diff = option_price - price
            vega_v = option.vega()

            # f(x) / f'(x)
            guess += diff / vega_v
        
        return guess    


    # according to the slides in W6:
    # These parameters - when fed back into the SVI expression above - will return a "variance" for a given strike.
    def SVI(self, a, b, p, m, sigma, KS):
        w = a + (b * ((p * (KS - m)) + np.sqrt(((KS - m)**(2)) + (sigma**(2)))))

        # As a reminder, volatility is expressed in terms of standard deviation.
        return np.sqrt(w)
    

    # the distances between SVI and implied volitilities
    def error(self, a, b, p, m, sigma, KS, implied_vol_vector):
        option = bs(self.S, self.K, self.r, self.q, self.sigma, self.T, self.option_type)
        SVI = option.SVI(a, b, p, m, sigma, KS)
        error = (SVI - implied_vol_vector)**2
        return np.mean(error)
    

    # constraints for optimization
    def condition_b(self, x):
        return x[1] # b ≥ 0
    
    def condition_rho(self, x):
        return 1 - abs(x[2]) # |ρ| < 1
    
    def condition_sigma(self, x):
        return x[4] # σ > 0
    

    # minimize the error function and get a, b, ρ, m, σ
    def mini(self, a, b, p, m, sigma, KS, implied_vol_vector):
    
        option = bs(self.S, self.K, self.r, self.q, self.sigma, self.T, self.option_type)
        obj_fun = lambda x: option.error(x[0], x[1], x[2], x[3] ,x[4], KS, implied_vol_vector)
        
        cons = [{'type':'ineq', 'fun':self.condition_b},
                {'type':'ineq', 'fun':self.condition_rho},
                {'type':'ineq', 'fun':self.condition_sigma}]
        para = [a, b, p, m, sigma]
        res = minimize(obj_fun, para, constraints = cons)

        return res.x