import numpy as np
import math
import matplotlib.pyplot as plt


gamma = float(input("Exhaust gases specific heats ratio(gamma):"))
M_exit = float(input("Exit gases Mach number:"))
N = int(input("Number of characteristics:"))
Rt = float(input("Throat Radius in cm:"))

def P_M_angle(mach, gam):  # P-M angle (Prandtl-Meyer angle)
        term1 = math.sqrt((gam - 1) / (gam + 1) * (mach ** 2 - 1))  # First term
        term2 = math.sqrt((gam + 1) / (gam - 1)) * math.atan(term1)  # Second Term
        term3 = math.atan(math.sqrt(mach ** 2 - 1))  # Third term
        nu = term2 - term3  # P-M angle in radians
        return math.degrees(nu)

#Get Mach number from gamma and v || Function
def get_mach(gam, nu):
        
        mach_guess = 1
        while True:
            nu_guess = P_M_angle(mach_guess, gam)
            if math.fabs(nu_guess - nu) < 0.2:
                break
            mach_guess += 0.005
        return mach_guess

#Initial variable entry
v_exit = P_M_angle(M_exit, gamma)
theta_max = v_exit/2
d_theta = theta_max/N

#Lists initialization
thetaList, vList, KplusList, KminusList, MList, muList, xList, yList = [], [], [], [], [], [], [], []

#Wall values list
xwallList, ywallList, thetawallList = [0], [Rt], [theta_max]




#Initial point characteristics: "0,n"

for i in range(N):

    #Values definiton
    theta = d_theta * (i + 1)
    v = theta
    M = get_mach(gamma, v)
    mu = math.degrees(math.sin(1/M))
    
    Kplus = theta - v
    Kminus = theta + v

    #Appending the values
    thetaList.append(theta)
    vList.append(v)
    MList.append(M)
    
    muList.append(mu)
    KplusList.append(Kplus)
    KminusList.append(Kminus)
    
    
k = 0

for i in range(1, N+1):

    for j in range(1, N-i+3):

        #First j point, on the symmetry line
        if j == 1:

            #Values definition
            theta = 0
            v = KminusList[i-1]
            M = get_mach(gamma, v)
            mu = math.degrees(math.sin(1/M))
            
            Kplus = theta - v
            Kminus = theta + v

            #Position definition
            slopeCminus = math.radians(0.5 * (thetaList[i-1] - muList[i-1] - mu))
            x = -1 * Rt/math.tan(slopeCminus)
            y = 0


            #Appending the values
            thetaList.append(theta)
            vList.append(v)
            MList.append(M)
            muList.append(mu)
            
            KplusList.append(Kplus)
            KminusList.append(Kminus)

            xList.append(x)
            yList.append(y)



        #Every other j point  
        elif 1 < j < N-i+2:

            #Values definition
            theta = 0.5 * (KminusList[j-1+k] + KplusList[-1])
            v = 0.5 * (KminusList[j-1+k] - KplusList[-1])
            M = get_mach(gamma, v)
            mu = math.degrees(math.sin(1/M))
            
            Kplus = theta - v
            Kminus = theta + v

            #Position definition
            slopeCminus = math.radians(0.5 * (thetaList[j-1+k] - muList[j-1+k] + theta - mu))
            slopeCplus = math.radians(0.5 * (thetaList[-1] + muList[-1] + theta + mu))

            x = (-xList[-1] * math.tan(slopeCplus) + yList[-1] - Rt)/(math.tan(slopeCminus) - math.tan(slopeCplus))
            y = x*math.tan(slopeCminus) + Rt


            #Appending the values
            thetaList.append(theta)
            vList.append(v)
            MList.append(M)
            muList.append(mu)
            
            KplusList.append(Kplus)
            KminusList.append(Kminus)

            xList.append(x)
            yList.append(y)


        #Last j point, wall point
        elif j == N-i+2:

            #Values definition
            theta = thetaList[-1]
            v = vList[-1]
            M = get_mach(gamma, v)
            mu = math.degrees(math.sin(1/M))
            
            Kplus = KplusList[-1]
            Kminus = KminusList[-1]

            #Setup linear equation solution with numpy
            A = np.array([
                [1, -math.tan(math.radians(theta + mu))],
                [1, -math.tan(math.radians(0.5*(theta + thetawallList[-1])))]
            ])

            a = -xList[-1] * math.tan(math.radians(theta + mu)) + yList[-1]
            b = -xwallList[-1] * math.tan(math.radians(0.5*(theta + thetawallList[-1]))) + ywallList[-1]
            
            B = np.array([a, b])

            #Solution
            solution = np.linalg.solve(A, B)
            x = solution[1]
            y = solution[0]


            #Appending the values
            thetaList.append(theta)
            vList.append(v)
            MList.append(M)
            muList.append(mu)
            
            KplusList.append(Kplus)
            KminusList.append(Kminus)

            #Wall values appending
            thetawallList.append(theta)
            xwallList.append(x)
            ywallList.append(y)

    k += 1


    
plt.figure(figsize=(7, 7))
plt.plot(xwallList, ywallList, linestyle='-', color='black', label='Nozzle Contour')
plt.title('Nozzle Contour')
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.show()




            
    
