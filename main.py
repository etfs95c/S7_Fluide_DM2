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
    return Beta

def beta_after_deviation(M_1, P_1, Theta, w):
    Beta = beta_with_M_Theta(M_1, Theta)[w]

    gamma = 1.4

    Mn_1 = M_1*np.sin(np.radians(Beta))
    Mn_2 = np.sqrt((1 + (gamma-1)/2) * Mn_1**2)/(gamma*Mn_1**2 - (gamma-1)/2)
    M_2 = (Mn_2/(np.sin(np.radians(Beta-Theta))))

    P_2 = P_1*(1 + (2*gamma)/(gamma+1)*(Mn_1**2 - 1))
    print("\nFor \u03B8 = ", Theta,"°, and M\u2081 =", M_1,"and P\u2081 =", M_1,"atm, we get\n \u03B2 = ", Beta,"°\tM\u2082 = ", M_2,"\tP\u2082 = ", P_2,"atm")
    return M_2, P_2


print("Hello World")

beta_after_deviation(100, 1, 20, 1)

print("Hello World")