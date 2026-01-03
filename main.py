import numpy as np
import math
import sys

def beta_with_M_Theta(M_1, Theta):
    Theta = math.radians(Theta)
    gamma = 1.4
    coef = [(gamma-1+(2/M_1**2))*math.tan(Theta), (-2+(2/M_1**2)), (gamma+1+(2/M_1**2))*math.tan(Theta)]
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

print("Hello World")

beta_with_M_Theta(100,20)
print("Hello World")