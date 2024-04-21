# Assignment

Your assignment, due in two weeks, is two-part. Part one is an
implementation of an implied volatility finder using both bisection and
Newton’s method. The user should be able to provide underlying price,
strike price, risk-free rate, tenor, option price, and call/put type. The code
should return the implied volatility measure.

Part two is an implementation of the Gatheral SVI skew fit algorithm.
Students must find or create a non-flat volatility skew - that is, a series of
implied volatility points - and write code that returns the five Gatheral SVI
parameters that represent a reasonable skew fit to the empirical data.
Results should resemble that which was presented in the lecture slides.