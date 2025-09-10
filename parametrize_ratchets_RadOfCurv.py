import numpy as np
from scipy.integrate import quad
from scipy.optimize import root

# In this program, R2 and R_curv correspond to the arclengths of the ratchet 
# in the xz and xy planes. 

#%% Compute a, b, c, from alpha, R2, R_curv, theta
def F(vars, alpha, R2, R_curv):
    a, b, c = vars

    # Eqn 1
    f1 = np.arctan(b/a) - alpha

    # Eqn 2
    f2 = c**2/a - R2

    # Eqn 3
    f3 = b**2/a - R_curv

    return np.array([f1, f2, f3])

# Example parameters
alpha = 60*np.pi/180
R2 = 0.4
R_curv = 0.75

# Initial guess

x0 = [1, 1, 1]

sol = root(F, x0, args=(alpha, R2, R_curv))

print("Success or Failure:", sol.success)
print("Solution (a, b, c):", sol.x[0], sol.x[1], sol.x[2])

a = sol.x[0]
b = sol.x[1]
c = sol.x[2]

#%% Compute alpha, R1, R2 for given values of a, c, theta

a = 40
b = 20
c = 45
theta = 40*np.pi/180

alpha = np.arctan(c/a)*180/np.pi

integrand1 = lambda z: np.sqrt(1 + (a**2/c**2 - 1) * np.sin(z)**2)
I1, _ = quad(integrand1, np.pi/2, 0)
S1 = -c * I1

S2 = a*theta

