import numpy as np

# reference Python for Finance p.295
# reference Hull ch.13
# reference Hull ch.21

# option payoff function
def payoff_function(option_type, S, K, 
                    alternative = 0):
    
    if option_type == "Call": 
        return np.maximum(S - K, alternative)
    
    elif option_type == "Put": 
        return np.maximum(K - S, alternative)
    
    else:
        raise ValueError("'Call' or 'Put'" )

# creat 2 empty lists to store price trees for Greeks calculation
for_greeks = []
STK = []

def binomial_tree(S, K, T, r, sigma, q, N, 
                  option_style, option_type):

    ######################## This is part one ########################

    # references Python for Finance p.295
    
    # dt, u, d, a, p, and discount factor
    dt = T / N
    u = np.exp(sigma * np.sqrt(dt))
    d = np.exp(-sigma * np.sqrt(dt)) 
    a = np.exp((r - q) * dt)
    p = (a - d) / (u - d)
    discount = np.exp(-r * dt)
    
    # ndarray object with gross upward movements
    _ = np.arange(N + 1)
    up = np.resize(_, (N + 1, N + 1))

    #  ndarray object with gross downward movements
    down = up.T * 2

    ST = S * np.exp(sigma * np.sqrt(dt) * (up - down))

    # append S tree to STK for Greeks calculation
    STK.append(ST)

    # get the last node and pass it to the recursive function
    last_node = payoff_function(option_type, ST[:,N], K)


    ######################## This is part two - Recursive Function ########################


    def binomial_recursive(payoff):
        
        if option_style == "European":

            pu = payoff[:-1]
            pd = payoff[1:]
            price = ((pu * p) + (pd * (1 - p))) * discount

            # store the first 3 nodes
            if price.shape[0] < 4:
                for_greeks.append(price)

            # return the price
            if price.shape[0] == 1:
                return price
            
        elif option_style == "American":
            
            pu = payoff[:-1]
            pd = payoff[1:]
            price = ((pu * p) + (pd * (1 - p))) * discount

            # generalize the col
            col = -(N - price.shape[0] + 2)


            stock_price = ST[:price.shape[0], col]

            # compare intrinsic value
            price = payoff_function(option_type, 
                                    S = stock_price, 
                                    K = K, alternative = price)

            # store the first 3 nodes
            if price.shape[0] < 4:
                for_greeks.append(price)
            
            # return the price
            if price.shape[0] == 1:
                return price            

        else:
            raise ValueError("'Call' or 'Put'" )
            
            
        return binomial_recursive(price)   
    
    # pass the last node to the recursive function
    output = binomial_recursive(last_node)

    return float(output)

######################## This is part three - Greeks Function ########################

# reference Hull ch.21

def price_and_greeks(S, K, T, r, sigma, q, N, 
                     option_style, option_type):

    # run binomial_tree to get price, C, and S
    price = binomial_tree(S, K, T, r, sigma, q, N, 
                          option_style, option_type)

    print(f"{option_style} {option_type} option")
    print(f"price: {price}")

    dt = T / N

    # indx C and S we need for Greeks
    f00 = for_greeks[2]
    f11 = for_greeks[1][0]
    f10 = for_greeks[1][1]
    f22 = for_greeks[0][0]
    f21 = for_greeks[0][1]
    f20 = for_greeks[0][2]

    S11 = STK[0][0][1]
    S10 = STK[0][1][1]
    S22 = STK[0][0][2]
    S21 = STK[0][1][2]
    S20 = STK[0][2][2]
    
    
    delta = (f11 - f10) / (S11 - S10)
    print(f"Delta: {delta}")

    gamma = (((f22 - f21) / (S22 - S21)) - ((f21 - f20) / (S21 - S20)))\
          / (0.5 * (S22 - S20))
    print(f"Gamma: {gamma}")

    theta = float(f21 - f00) / (2 * dt)
    print(f"Theta: {theta}")

    # set delta sigma as 0.1 sigma
    ds = sigma / 10
    # C(sigma + ds)
    price_a =  binomial_tree(S, K, T, r, sigma + ds , q, 
                             N, option_style, option_type)
    # C(sigma - ds)
    price_b =  binomial_tree(S, K, T, r, sigma - ds , q, 
                             N, option_style, option_type)
    vega = (price_a - price_b) / (2 * ds)
    print(f"Vega: {vega}")

    #set delta r as 0.1 r
    dr = r / 10
    # C(r + dr)
    price_c =  binomial_tree(S, K, T, r + dr, sigma , q, 
                             N, option_style, option_type)
    #C(r - dr)
    price_d =  binomial_tree(S, K, T, r - dr, sigma , q, 
                             N, option_style, option_type)
    rho = (price_c - price_d) / (2 * dr)
    print(f"Rho: {rho}")

# Hull Ch.21.1 example
S = 50
K = 50
r = 0.1
sigma = 0.4
T = 5 / 12
q = 0
N = 5
option_style = "American"
option_type = "Put"
price_and_greeks(S, K, T, r, sigma, q, N, option_style, option_type)