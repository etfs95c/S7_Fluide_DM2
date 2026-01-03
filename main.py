import numpy as np
import math
import sys

def beta_with_M_Theta(M_1, Theta):
    Theta_rad = np.radians(Theta)
    gamma = 1.4
    coef = [(gamma-1+(2/M_1**2))*np.tan(Theta_rad), (-2+(2/M_1**2)), (gamma+1+(2/M_1**2))*np.tan(Theta_rad)]
    print("\nHello1 \n")

    roots = np.roots(coef)

    if roots.dtype == complex:
        print("One or more root are complex, no physical solution")
        sys.exit(1)
    elif len(roots) <=1:
        print("No physical solution")
        sys.exit(2)
    else:
        Beta = np.degrees(np.arctan(roots))
    
    print("The value of \u03B2 for a strong shock is :\n\t\u03B2 = ", Beta[0] ,"\nThe value of \u03B2 for a weak shock is :\n\t\u03B2 = ", Beta[1])
    print(roots)

def beta_after_deviation(M_1, P_1, theta):
    Beta_s = beta_with_M_Theta[0]
    Beta_w = beta_with_M_Theta[1]

    gamma = 1.4

    Mn_1_s = M_1*np.sin(np.radians(Beta_s))
    Mn_2_s = np.sqrt((1 + (gamma-1)/2) * Mn_1_s**2)/(gamma*Mn_1_s**2 - (gamma-1)/2)

print("Hello World")

beta_with_M_Theta(100,20)
print("Hello World")