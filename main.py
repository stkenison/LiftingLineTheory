import numpy as np; import os;import matplotlib.pyplot as plt; import json
os.system('cls'); plt.rcParams['font.family'] = 'Times New Roman'

#get inputs from .json file
with open('input.json') as f: data = json.load(f)

#Load base data from output file
R_A = data['wing']['planform']['aspect_ratio'] #Load aspect ratio from file
R_T = data['wing']['planform']['taper_ratio'] #Load taper ratio from file
N = data['wing']['nodes_per_semispan']*2+1 #Load desired semispan nodes and calculate total nodes
C_La_2D = data['wing']['airfoil_lift_slope'] #Load 2D airfoil lift slope

#Cosine clustering to form airfoil
theta = np.linspace(0,np.pi,N) #calculate theta values based off number of nodes
z_b = -0.5*np.cos(theta)
c_b = 2*(1-(1-R_T)*np.abs(np.cos(theta)))/(R_A*(1+R_T)) #chord length normalized by span

#plot airfoil for verification
plt.figure('Planform Plot'); plt.plot(z_b,c_b/4, color = 'blue')
plt.plot(z_b,-3*c_b/4, color = 'blue')
plt.xlabel("Spanwise Position (z/b)")
plt.ylabel("Chord Position (c/b)")
plt.plot([z_b[0],z_b[0]],[-3*c_b[0]/4,c_b[0]/4],color = 'blue'); plt.plot([z_b[-1],z_b[-1]],[-3*c_b[-1]/4,c_b[-1]/4],color = 'blue')
for i in range(N):
    plt.plot([z_b[i],z_b[i]],[-3*c_b[i]/4,c_b[i]/4],color = 'red', linestyle = '--')
plt.axis("equal")
plt.grid(True)
plt.show()

washout_distribution_type = data['wing']['washout']['distribution'] #get washout distribution type
if washout_distribution_type == 'optimum':
    omega = 1-np.sin(theta)/(c_b/c_b[N//2]) #define optimum washout distribution
elif washout_distribution_type == 'linear':
    omega = 1-np.sin(theta)
else: 
    omega = np.zeros(N)
washout = np.radians(data['wing']['washout']['amount[deg]']) #Load washout amount from input file

C = np.zeros((N,N)) #initialize C array
for i in range(N):
    for j in range(N):
        n = j+1
        if i==0: 
            C[i,j] = n**2
        elif i==N-1:
            C[i,j] = ((-1)**(n+1))*(n**2)
        else:
            C[i,j] = (4/C_La_2D/c_b[i]+n/np.sin(theta[i]))*np.sin(n*theta[i])
C_inv = np.linalg.inv(C) #POSSIBLE ISSUES
a = np.matmul(C_inv,np.ones(N))
b = np.matmul(C_inv,omega)

#calculate kappa_L value
kappa_L = (1-(1+R_A*np.pi/C_La_2D)*a[0])/((1+R_A*np.pi/C_La_2D)*a[0])

#Calculate C_La_3D value
C_La_3D = np.pi*R_A*a[0]

#calculate epsilon_omega value
epsilon_omega = b[0]/a[0]
print(epsilon_omega) ##COULD HAVE ISSUES

#calculate kappa_D value
kappa_D = 0 #initialize kappa_D
for i in range(1,N): #sum kappa_D values
    kappa_D += (i+1)*(a[i]/a[0])**2 #calculate part of sum

#calculate kappa_DL value
kappa_DL = 0 #initialize kappa_DL
for i in range(1,N): #sum kappa_DL values
    kappa_DL += 2*b[0]/a[0]*(i+1)*a[i]/a[0]*(b[i]/b[0]-a[i]/a[0]) #calculate part of sum

#calculate kappa_Domega value
kappa_Domega = 0 #initialize kappa_Domega
for i in range(1,N): #sum kappa_Domega values
    kappa_Domega += (b[0]/a[0])**2*(i+1)*(b[i]/b[0]-a[i]/a[0])**2 #calculate part of sum

#calculate kappa_D0 value
kappa_D0 = kappa_D-(kappa_DL**2/(4*kappa_Domega))

with open('output.txt', 'w') as f:
    f.write('C = \n\n'); np.savetxt(f, C, delimiter=', ', fmt='%.5f')
    f.write('\n\n\nC_inv = \n\n'); np.savetxt(f, C_inv, delimiter=', ', fmt='%.5f')
    f.write('\n\n\na_n = \n\n'); np.savetxt(f, a, delimiter=', ', fmt='%.5f')
    f.write('\n\n\nb_n = \n\n'); np.savetxt(f, b, delimiter=', ', fmt='%.5f')




'''
#initialize variables
geometry = np.loadtxt(data['geometry'],dtype = float) #Load airfoil geometry from file
alpha = np.array(data['alpha[deg]'])/180*np.pi #define desired angles of attack
if np.isscalar(alpha): alpha = np.array([alpha])  # Convert scalar to an array if singular value
v_inf = data['freestream_velocity'] #define freestream velocity
n = geometry.shape[0] #count number of points in airfoil geometry file
n_alpha = len(alpha) #count number of unique values in angle of attack list

#initialize variables to describe geometry
CP = 0.5*(geometry[1:]+geometry[:-1])#define control points
dx = geometry[1:,0]-geometry[:-1,0] #define dx values
dy = geometry[1:,1]-geometry[:-1,1] #define dy values
l = np.sqrt(dx**2+dy**2) #define lengths of panels

#solve for A array
A = np.zeros((n,n)) #initialize A array
for i in range(n-1):
    for j in range(n-1):
        xi = (dx[j]*(CP[i,0]-geometry[j,0])+dy[j]*(CP[i,1]-geometry[j,1]))/l[j] #define xi term
        eta = (-dy[j]*(CP[i,0]-geometry[j,0])+dx[j]*(CP[i,1]-geometry[j,1]))/l[j] #define eta term
        
        phi = np.arctan2(eta*l[j],eta**2+xi**2-xi*l[j]) #define phi term
        psi = 0.5*np.log((xi**2+eta**2)/((xi-l[j])**2+eta**2)) #define psi term

        P1 = np.array([[dx[j],-dy[j]],[dy[j],dx[j]]]) #calculate first part of P array
        P2 = np.array([[(l[j]-xi)*phi+eta*psi,xi*phi-eta*psi],
            [eta*phi-(l[j]-xi)*psi-l[j],-eta*phi-xi*psi+l[j]]]) #calculate second part of P array
        P = (1/(2*np.pi*l[j]**2))*np.matmul(P1,P2) #create P array by matrix multiplying parts

        A[i,j] += dx[i]*P[1,0]/l[i]-dy[i]*P[0,0]/l[i] #calculat current cell of A array
        A[i,j+1] += dx[i]*P[1,1]/l[i]-dy[i]*P[0,1]/l[i] #start calculating next cell of A array

A[n-1,0] = 1.0 #define edge of A array based on kutta condition
A[n-1,n-1] = 1.0 #define corner of A array based on kutta condition

#use B array to solve for gamma for all angles of attack
B = np.zeros((n,n_alpha)) #initialize B column(s)
gamma = np.zeros((n,n_alpha))  #initialize gamma column(s)
for i in range(n_alpha):
    B[:-1,i] = v_inf*((geometry[1:,1]-geometry[:-1,1])*np.cos(alpha[i])-
        (geometry[1:,0]-geometry[:-1,0])*np.sin(alpha[i]))/[l[:]] #create B array for current angles of attack
    gamma[:,i] = np.linalg.solve(A,B[:,i]) #solve for gamma asoociated with current angle of attack

#find coeffecient of lift
C_L = np.zeros(n_alpha) #initialize lift coefficient array
for i in range(n_alpha):
    C_L[i] = np.sum(l[:]*(gamma[:-1,i]+gamma[1:,i])/v_inf) #find value of current iteration of lift coefficient array

#find leading edge moment coefficient array
C_mLE = np.zeros(n_alpha) #initialize leading edge moment coefficient array
for i in range(n_alpha):
    C_mLE[i] = -1/3*np.sum(l[:]/v_inf*((2*geometry[:-1,0]*gamma[:-1,i]+geometry[1:,0]*gamma[:-1,i]+geometry[:-1,0]*gamma[1:,i]+2*geometry[1:,0]*gamma[1:,i])*np.cos(alpha[i])
    +(2*geometry[:-1,1]*gamma[:-1,i]+geometry[1:,1]*gamma[:-1,i]+geometry[:-1,1]*gamma[1:,i]+2*geometry[1:,1]*gamma[1:,i])*np.sin(alpha[i]))) #find value of current iteration of lift coefficient array

# Find quarter chord moment coefficient array
C_mQC = np.zeros(n_alpha)  #initialize quarter chord moment coefficient array
QC_x = geometry[:,0]-0.25 #Compute the moment arm relative to the quarter chord
for i in range(n_alpha):
    C_mQC[i] = -1/3 * np.sum(
        l[:]/v_inf*(
            (2*QC_x[:-1]*gamma[:-1, i]+QC_x[:-1]*gamma[1:,i]+ 
            QC_x[1:]*gamma[:-1,i]+2*QC_x[1:]*gamma[1:,i])*np.cos(alpha[i])+
            (2*geometry[:-1,1]*gamma[:-1,i]+geometry[1:,1]*gamma[:-1,i]+ 
            geometry[:-1,1]*gamma[1:,i]+2*geometry[1:,1]*gamma[1:,i])*np.sin(alpha[i])
        )
    )  # find value of current iteration of quarter chord moment coefficient array

#print out calculated coefficients for verification
print(f"a: {alpha*180/np.pi}\n")
print(f"C_L: {C_L}\n")
print(f"C_mLE: {C_mLE}\n")
print(f"C_mQC: {C_mQC}\n")

#Plot airfoil
plt.figure(1)
plt.plot(geometry[:,0],geometry[:,1])
plt.axis("equal")
plt.title('Airfoil Geometry')

#Plot coefficient of lift
plt.figure(2)
plt.plot(alpha*180/np.pi,C_L)
plt.title('Lift Coefficient ($C_{L}$)')
plt.xlabel("Angle ($^{o}$)");plt.ylabel("$C_{L}$",rotation='horizontal')

#Plot leading edge moment coefficient
plt.figure(3)
plt.plot(alpha*180/np.pi,C_mLE)
plt.title('Coefficient of Leading Edge Moment ($C_{mLE}$)')
plt.xlabel("Angle ($^{o}$)");plt.ylabel("$C_{mLE}$",rotation='horizontal')

#Plot quarter-chord moment coefficient
plt.figure(4)
plt.plot(alpha*180/np.pi,C_mQC)
plt.title('Coefficient of Quarter-Chord Moment ($C_{mc/4}$)')
plt.xlabel("Angle ($^{o}$)");plt.ylabel("$C_{mc/4}$",rotation='horizontal')
plt.show()
'''