import numpy as np
from scipy.integrate import quad
from scipy.optimize import root


#%% Compute alpha, R1, R2 from a, b, c, theta

a = 2
b = 4.23
c = 2.58

theta = 20*np.pi/180

alpha = np.arctan(c/a)*180/np.pi

integrand1 = lambda z: - np.sqrt(1 + (a**2/c**2 - 1) * np.sin(z)**2)
R1, _ = quad(integrand1, np.pi/2, 0)


integrand2 = lambda z: -np.sqrt(1 + (a**2/b**2 - 1) * np.sin(z)**2)
R2, _ = quad(integrand2, np.pi/2 + theta/2, np.pi/2 - theta/2)


print(f"alpha = {alpha}, R1 = {R1}, R2 = {R2}")

#%% Compute a, b, c, from alpha, R1, R2, theta
def F(vars, alpha, R1, R2, theta):
    a, b, c = vars

    # Eqn 1
    f1 = np.arctan(c/a) - alpha

    # Eqn 2
    integrand1 = lambda z: np.sqrt(1 + (a**2/c**2 - 1) * np.sin(z)**2)
    I1, _ = quad(integrand1, np.pi/2, 0)
    f2 = -c * I1 - R1

    # Eqn 3
    integrand2 = lambda z: np.sqrt(1 + (a**2/b**2 - 1) * np.sin(z)**2)
    I2, _ = quad(integrand2, np.pi/2 + theta/2, np.pi/2 - theta/2)
    f3 = -b * I2 - R2

    return np.array([f1, f2, f3])

# Example parameters
# alpha = 60.8*np.pi/180
# R1 = 1.248
# R2 = 0.123
theta = 20*np.pi/180

# Initial guess
result = "False"
while result == "False":
    
    x0 = np.random.uniform(0.1, 10, 3)
    x0 = [2,4,3]

    sol = root(F, x0, args=(alpha, R1, R2, theta))

    print("Success or Failure:", sol.success)
    print("Solution (a, b, c):", sol.x)

a = sol.x[0]
b = sol.x[1]
c = sol.x[2]
