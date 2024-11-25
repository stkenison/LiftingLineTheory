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

#define washout distribution
washout_distribution_type = data['wing']['washout']['distribution'] #get washout distribution type
if washout_distribution_type == 'optimum':
    omega = 1-np.sin(theta)/(c_b/c_b[N//2]) #define optimum washout distribution
elif washout_distribution_type == 'linear':
    omega = np.abs(np.cos(theta)) #define linear washout distribution, eqn 6.30
else: 
    omega = np.zeros(N)

#define washout based off of inputs
if type(data['wing']['washout']['amount[deg]']) == int or type(data['wing']['washout']['amount[deg]']) == float:
    washout_amount = np.radians(data['wing']['washout']['amount[deg]']) #Load washout amount from input file
else:
    washout_amount = 2*(R_T+1)*data['wing']['washout']['CL_design']/(np.pi*C_La_2D) #define washout based on optimum washout distribution for wings with linear taper
print('washout = ', washout_amount)

#create C and C inverse arrays
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
C_inv = np.linalg.inv(C)

#solve for a and b matrixes used for coefficient solutions
a = np.matmul(C_inv,np.ones(N))
b = np.matmul(C_inv,omega)

#define theta_f for aileron distribution
theta_f = np.zeros(N)
begin_z_b = data['wing']['aileron']['begin[z/b]']
end_z_b = data['wing']['aileron']['end[z/b]']
begin_cf_c = data['wing']['aileron']['begin[cf/c]']
end_cf_c = data['wing']['aileron']['end[cf/c]']
for i in range(N):
    if (z_b[i] > begin_z_b and z_b[i] < end_z_b):
        cf_c = begin_cf_c+(abs(z_b[i])-begin_z_b)*(end_cf_c-begin_cf_c)/(end_z_b-begin_z_b)
        theta_f[i] = np.arccos(2*cf_c-1)
    elif (-z_b[i] > begin_z_b and -z_b[i] < end_z_b):
        cf_c = begin_cf_c+(abs(z_b[i])-begin_z_b)*(end_cf_c-begin_cf_c)/(end_z_b-begin_z_b)
        theta_f[i] = np.arccos(2*cf_c-1)

#define aileron distribution
chi = np.zeros(N) #initialize aileron distribution
for i in range(N):
    if (z_b[i] > begin_z_b and z_b[i] < end_z_b):
        chi[i] = -data['wing']['aileron']['hinge_efficiency']*(1-(theta_f[i]-np.sin(theta_f[i]))/(np.pi))
    elif (-z_b[i] > begin_z_b and -z_b[i] < end_z_b):
        chi[i] = data['wing']['aileron']['hinge_efficiency']*(1-(theta_f[i]-np.sin(theta_f[i]))/(np.pi)) 

#calculate c_n matrix
c = np.matmul(C_inv,chi)

#calculate d_n matrix
d = np.matmul(C_inv,np.cos(theta))

#calculate kappa_L value
kappa_L = (1-(1+R_A*np.pi/C_La_2D)*a[0])/((1+R_A*np.pi/C_La_2D)*a[0])

#Calculate C_La_3D value
C_La_3D = np.pi*R_A*a[0]

#calculate epsilon_omega value
epsilon_omega = b[0]/a[0] ##COULD HAVE ISSUES

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

#calculate C_lp value
C_lda = -np.pi*R_A*c[1]/4

#calculate C_lp value 
C_lp = -np.pi*R_A*d[1]/4

#calculate coefficient of lift
if type(data['condition']['alpha_root[deg]']) == int or type(data['condition']['alpha_root[deg]']) == float:
    alpha_root = np.radians(data['condition']['alpha_root[deg]'])
    C_L = C_La_3D*(alpha_root-epsilon_omega*washout_amount)
else:
    C_L = data['wing']['condition']['CL']
    alpha_root = np.radians(1) #calculate root angle of attack

#calculate induced drag coefficient
C_Di = 1

#calculate rolling moment coefficient
C_l = 1

#calculate yawing moment coefficient
C_n = 1

#calculate steady dimenionless roll rate
p_steady = 1



#write output matrixes to external file
with open('output.txt', 'w') as f:
    f.write('C = \n\n'); np.savetxt(f, C, delimiter=', ', fmt='%.5f')
    f.write('\n\n\nC_inv = \n\n'); np.savetxt(f, C_inv, delimiter=', ', fmt='%.5f')
    f.write('\n\n\na_n = \n\n'); np.savetxt(f, a, delimiter=', ', fmt='%.7f')
    f.write('\n\n\nb_n = \n\n'); np.savetxt(f, b, delimiter=', ', fmt='%.7f')
    f.write('\n\n\nc_n = \n\n'); np.savetxt(f, c, delimiter=', ', fmt='%.7f')
    f.write('\n\n\nd_n = \n\n'); np.savetxt(f, d, delimiter=', ', fmt='%.7f')

#print calculated constant values to terminal
print("K_L =",kappa_L)
print("C_L,\u03B1 =",C_La_3D)
print("\u03B5_\u03A9 =", epsilon_omega,"COULD HAVE ISSUES") 
print("K_D =",kappa_D)
print("K_DL =",kappa_DL)
print("K_D\u03A9 = =",kappa_Domega)
print("K_Do =",kappa_D0)
print("C_l,\u03B4\u03B1 =",C_lda)
print("C_l,\u0070\u0304 =",C_lp,"\n")

#print calculated perfomance values to terminal
print("C_L =",C_L)
print("C_Di",C_Di)
print("C_l =",C_l)
print("C_n =",C_n)
print("\u0070\u0304 =",p_steady,"\n")

if data['view']['planform']: #plot planform if desired by user
    plt.figure('Planform Plot'); plt.plot(z_b,c_b/4, color = 'blue')
    plt.plot(z_b,-3*c_b/4, color = 'blue')
    plt.xlabel("Spanwise Position (z/b)")
    plt.ylabel("Chord Position (c/b)")
    plt.plot([z_b[0],z_b[0]],[-3*c_b[0]/4,c_b[0]/4],color = 'blue'); plt.plot([z_b[-1],z_b[-1]],[-3*c_b[-1]/4,c_b[-1]/4],color = 'blue')
    for i in range(N):
        plt.plot([z_b[i],z_b[i]],[-3*c_b[i]/4,c_b[i]/4],color = 'red', linestyle = '--')
    plt.axis("equal")
    plt.grid(True)
    right_aileron = np.array([[begin_z_b,end_z_b,end_z_b,begin_z_b,begin_z_b],[-1/R_A*3/4,-1/R_A*3/4,-1/R_A*3/4*(1-end_cf_c),-1/R_A*3/4*(1-begin_cf_c),-1/R_A*3/4]])
    plt.plot(right_aileron[0,:],right_aileron[1,:],color = 'green')
    plt.plot(-right_aileron[0,:],right_aileron[1,:],color = 'green') #flip aileron coordinate array for left side

if data['view']['washout_distribution']: #plot washout distribution if desired by user
    plt.figure('Washout Distribution')
    plt.plot(z_b,omega, color = 'blue')
    plt.grid(True)
    plt.ylabel("Washout Distribution (\u03C9(z))")
    plt.xlabel("Spanwise Position (z/b)")

if data['view']['aileron_distribution']: #plot aileron distribution if desired by user
    plt.figure('Aileron Distribution')
    plt.plot(z_b,chi, color = 'blue')
    plt.grid(True)
    plt.ylabel("Aileron Position (\u03A7(z))")
    plt.xlabel("Spanwise Position (z/b)")

if data['view']['CL_hat_distributions']: #plot CL hat distributions if desired by user
    plt.figure('CL_hat_distributions')

if data['view']['CL_tilde_distributions']: #plot CL tilde distributions if desired by user
    plt.figure('CL_tilde_distributions')

plt.show()
