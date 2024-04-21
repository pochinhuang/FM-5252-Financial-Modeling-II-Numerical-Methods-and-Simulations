# Week 8 Lecture
# Random in C#

The first transformation algorithm we’ll call **sum twelve**. This algorithm
works in the following way.

1. Generate twelve uniform random numbers on the unit interval

2. Add these twelve values together

3. Subtract six from the total

Although perhaps unintuitive, the output is a Gaussian random variable.
In practice, this is a poor algorithm due to computational inefficiency and
statistical gremlins.


---


The second transformation algorithm is known as the **Box-Muller**
**transform**. It is implemented as follows.

1. Generate two uniform random values on the unit interval

2. Input them to the following expressions

3. The resulting output is a pair of Gaussian random variables

$z_1 = \sqrt{- 2\ln(x_1)} \cos(2\pi x_2)$
$z_2 = \sqrt{- 2\ln(x_1)} \sin(2\pi x_2)$


---


The final transform is known as **polar rejection**. It works as follows.

1. Generate two uniform random values, $x_1$, $x_2$
2. w = $x_{1}^{2} + x_{2}^{2}$
3. If w >1, repeat first two steps, otherwise...
4. c = $\sqrt{- 2\ln(w) / w}$
5. $z_1 = cx_1$
6. $z_2 = cx_2$