import numpy as np
import math
import sys
import matplotlib.pyplot as mplt
from matplotlib.ticker import MultipleLocator

def beta_with_M_Theta(M_1, Theta):
    Theta_rad = np.radians(Theta)
    gamma = 1.4
    coef = [(gamma-1+(2/M_1**2))*np.tan(Theta_rad), (-2+(2/M_1**2)), (gamma+1+(2/M_1**2))*np.tan(Theta_rad), 2/(M_1**2)]

    roots = np.roots(coef)

    if roots.dtype == complex:
        print("One or more root are complex, no physical solution")
        print(roots)
        sys.exit(1)
    elif len(roots) <=1:
        print("No physical solution")
        print(roots)
        sys.exit(2)
    else:
        Beta = np.degrees(np.arctan(roots))
    
    #print("The value of \u03B2 for a strong shock is :\n\t\u03B2 = ", round(Beta[0],3) ,"\nThe value of \u03B2 for a weak shock is :\n\t\u03B2 = ", round(Beta[1],3))
    #print(roots)
    #Desactivated after Q3.

    return Beta

def beta_after_deviation(M_1, P_1, Theta, w): #w=1 for weak shock, w=0 for strong shock
    Beta = round(beta_with_M_Theta(M_1, Theta)[w],3)

    gamma = 1.4

    Mn_1 = M_1*np.sin(np.radians(Beta))
    Mn_2 = np.sqrt((1 + (gamma-1)/2 * Mn_1**2)/(gamma*Mn_1**2 - (gamma-1)/2))
    M_2 = round((Mn_2/(np.sin(np.radians(Beta-Theta)))),3)

    P_2 = round(P_1*(1 + (2*gamma)/(gamma+1)*(Mn_1**2 - 1)),3)

    #print("\nFor \u03B8 = ", Theta,"°, and M\u2081 =", M_1,"and P\u2081 =", M_1,"atm, we get\n \u03B2 = ", Beta,"°\tM\u2082 = ", M_2,"\tP\u2082 = ", P_2,"atm")
    #Desactivated after Q3.

    return[M_2, P_2, Beta]

Flow_1 = [3,1]
print("M\u2081 = ", Flow_1[0],"\tP\u2081 = ", Flow_1[1]," atm\n")

temp = beta_after_deviation(Flow_1[0], Flow_1[1], 20, 1)
Flow_2 = [temp[0],temp[1]]
Beta_2 = temp[2]
print("M\u2082 = ", Flow_2[0],"\tP\u2082 = ", Flow_2[1]," atm\t\u03B2\u2082 = ", Beta_2,"°\n")

temp = beta_after_deviation(Flow_1[0], Flow_1[1], 15, 1)
Flow_3 = [temp[0],temp[1]]
Beta_3 = temp[2]
print("M\u2083 = ", Flow_3[0],"\tP\u2083 = ", Flow_3[1]," atm\t\u03B2\u2083 = ", Beta_3,"°\n")

del temp

#Diagram State 1
n = 100
theta_range = np.linspace(0, 34.073,n)
P_vect_weak = []
P_vect_strong = []
for i in theta_range:
    P_vect_weak.append(beta_after_deviation(Flow_1[0], Flow_1[1], i, 1)[1])
    P_vect_strong.append(beta_after_deviation(Flow_1[0], Flow_1[1], i, 0)[1])
P_vect_strong[0]=P_vect_strong[1]

theta_range_1 = np.concatenate([-np.array(theta_range)[::-1], np.array(theta_range)])
P_vect_weak_1 = np.concatenate([np.array(P_vect_weak)[::-1], np.array(P_vect_weak)])
P_vect_strong_1 = np.concatenate([np.array(P_vect_strong)[::-1], np.array(P_vect_strong)])

mplt.figure(figsize=(10, 6))

mplt.plot(theta_range_1, P_vect_weak_1, label='Weak Shock Solution', color='deepskyblue', linestyle='-', linewidth=2)
mplt.plot(theta_range_1, P_vect_strong_1, label='Strong Shock Solution', color='limegreen', linestyle='-', linewidth=2)

mplt.title('Pressure Deflection Diagram State 1', fontsize=18)
mplt.xlabel('Deflection Angle θ (°)', fontsize=14)
mplt.ylabel('Static Pressure P (atm)', fontsize=14)
mplt.gca().xaxis.set_major_locator(MultipleLocator(10))  
mplt.gca().yaxis.set_major_locator(MultipleLocator(1))   

mplt.grid(True, linestyle='--', color='gray', alpha=0.3)
mplt.legend(loc='best', fontsize=12)
mplt.gcf().set_facecolor('lightgray')
mplt.tight_layout()
mplt.show()

#Diagram State 2
theta_range=np.linspace(0, 22.869,n)       
P_vect_weak=[]
P_vect_strong=[]
for i in theta_range:
    P_vect_weak.append(beta_after_deviation(Flow_2[0], Flow_2[1], i, 1)[1])
    P_vect_strong.append(beta_after_deviation(Flow_2[0], Flow_2[1], i, 0)[1])
P_vect_strong[0]=P_vect_strong[1]


theta_range_2=np.concatenate([-np.array(theta_range)[::-1], np.array(theta_range)])
P_vect_weak_2=np.concatenate([np.array(P_vect_weak)[::-1], np.array(P_vect_weak)])
P_vect_strong_2=np.concatenate([np.array(P_vect_strong)[::-1], np.array(P_vect_strong)])


mplt.figure(figsize=(10, 6))

mplt.plot(theta_range_2, P_vect_weak_2, label='Weak Shock Solution', color='deepskyblue', linestyle='-', linewidth=2)
mplt.plot(theta_range_2, P_vect_strong_2, label='Strong Shock Solution', color='limegreen', linestyle='-', linewidth=2)

mplt.title('Pressure Deflection Diagram State 2', fontsize=18)
mplt.xlabel('Deflection Angle θ (°)', fontsize=14)
mplt.ylabel('Static Pressure P (atm)', fontsize=14)
mplt.gca().xaxis.set_major_locator(MultipleLocator(5))  
mplt.gca().yaxis.set_major_locator(MultipleLocator(1))   

mplt.grid(True, linestyle='--', color='gray', alpha=0.3)
mplt.legend(loc='best', fontsize=12)
mplt.gcf().set_facecolor('lightgray')
mplt.tight_layout()
mplt.show()


#Diagram State 3
theta_range=np.linspace(0, 26.8605, n)       
P_vect_weak=[]
P_vect_strong=[]
for i in theta_range:
    P_vect_weak.append(beta_after_deviation(Flow_3[0], Flow_3[1], i, 1)[1])
    P_vect_strong.append(beta_after_deviation(Flow_3[0], Flow_3[1], i, 0)[1])
P_vect_strong[0]=P_vect_strong[1]


theta_range_3=np.concatenate([-np.array(theta_range)[::-1], np.array(theta_range)])
P_vect_weak_3=np.concatenate([np.array(P_vect_weak)[::-1], np.array(P_vect_weak)])
P_vect_strong_3=np.concatenate([np.array(P_vect_strong)[::-1], np.array(P_vect_strong)])

mplt.figure(figsize=(10, 6))

mplt.plot(theta_range_3, P_vect_weak_3, label='Weak Shock Solution', color='deepskyblue', linestyle='-', linewidth=2)
mplt.plot(theta_range_3, P_vect_strong_3, label='Strong Shock Solution', color='limegreen', linestyle='-', linewidth=2)

mplt.title('Pressure Deflection Diagram State 3', fontsize=18)
mplt.xlabel('Deflection Angle θ (°)', fontsize=14)
mplt.ylabel('Static Pressure P (atm)', fontsize=14)
mplt.gca().xaxis.set_major_locator(MultipleLocator(5))  
mplt.gca().yaxis.set_major_locator(MultipleLocator(1))   

mplt.grid(True, linestyle='--', color='gray', alpha=0.3)
mplt.legend(loc='best', fontsize=12)
mplt.gcf().set_facecolor('lightgray')
mplt.tight_layout()
mplt.show()