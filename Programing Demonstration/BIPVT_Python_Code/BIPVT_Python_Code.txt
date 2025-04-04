import math

#inputs%%%%%%%%%%%
loop=2
forced=1
hours=10
nu=0.000015
ro=1.200
mu=1.51*(10**-5)
v_amb=2
g=9.8
A_pv=1
b=1
C_p=1005
p_2=1.184
k_c=0.036
k_g=1
k_i=0.035
k_tg=0.033
k_air=0.025
L=1
L_c=300*(10**(-6))
L_g=0.003
L_i=0.05
L_tg=0.0005
alpha_c=0.8
beta_c=0.9
alpha_tg=0.5
theta=math.radians(0)
etha_c = [0.2] * hours
etha_th=0.5
t_g=0.95
pr=0.72
D_h=9/100
ductdepth=0.06
A_2=ductdepth*b
I=[300, 400, 550, 650, 658, 600, 560, 420, 250, 110]
T_r=[32, 33, 35 ,36, 42 ,44, 42, 41, 41, 39]
T_r=list(map(lambda x:x+273 ,T_r)) # C to K

T_amb=[29, 30, 32, 35, 38, 40, 41, 41, 40, 39]
T_amb=list(map(lambda x:x+273 ,T_amb)) # C to K
v=[2.67, 1.98, 1.42, 1.87, 1.73, 1.67, 1.77, 1.7, 1.83, 1.77]

#Backsurface, Cell Temperature, Outlet air Temperature Insulation Temperature Initial Guesses%%%%%%%
T_bs=[38, 43, 47, 51, 57, 56, 55, 49, 46, 42]
T_bs=list(map(lambda x:x+273 ,T_bs)) # C to K

T_c=[38, 41, 48, 50, 55, 55, 53, 51, 47, 42]
T_c=list(map(lambda x:x+273 ,T_c)) # C to K

T_2=[33, 36, 39, 42, 46, 46, 46, 45, 43, 40]
T_2=list(map(lambda x:x+273 ,T_2)) # C to K

T_i = [273 + 30] * hours

T_i=list(map(lambda x:x+273 ,T_i)) # C to K

###defining some variables
Ra=[0]*hours
Gr=[0]*hours
f=[0]*hours
Raa=[0]*hours
Re=[0]*hours
h_air=[0]*hours
gg=[107000000]*hours
u_tg=(L_tg/k_tg)**(-1)
Ra_2=[0]*hours
h_r=[0]*hours
alpha=[0]*hours
h_rad=[0]*hours
u_b=[0]*hours
u_T=[0]*hours
u_tair=[0]*hours
u_bb=[0]*hours
u_L=[0]*hours
mdot=[0]*hours
zzz=[0]*hours
T_airout=[0]*hours
T_bs=[0]*hours
T_c=[0]*hours
qdot=[0]*hours
electrical11=[0]*hours
thermal=[0]*hours
tinsulation=[0]*hours
electrical1=[0]*hours
#######

#Calculation%%%%%%%%%%%%

for w in range(loop):

    #(Ra)/((H)^3(deltaT))%%
    #Ra, Gr, f%%%%%%%%%%%%%%%%%%

    for i in range(10):
        Ra[i]=abs((T_bs[i]-T_i[i]))*gg[i]*ductdepth**3
        Gr[i]=Ra[i]/pr
        f[i]=1.906*(Gr[i]/pr)**(1/12)

            #h_air(Duct Air)%%%%%%%%%%%%%%%%%%%%%%%
        
        if forced==0:
            Raa[i]=abs((T_bs[i]-T_i[i]))*(gg[i]*nu*math.sin(theta))*ro*(L**3)
        
            #Natural Convection%%%%%%
            
            h_air[i]=(k_air/D_h)*0.0965*(Raa[i]**0.29)

        else:
            
            Re[i]=(p_2*v[i]*D_h)/(mu)

        if Re[i]>2300:
            #Forced Convection(Turbulent)%%%%%%    
            h_air[i]=(k_air/D_h)*0.023*(Re[i]**0.8)*(pr**0.4)
        else:
            #Forced Convection(Laminar)%%%%%%    
            h_air[i]= (k_air/D_h)*2.98    #.*ones(1,10)

    #h_r(Room Side Natural Convection)%%%%%%%%%%%%%%%%%%%
    gggg=1+((0.492/pr)**(9/16))
    for i in range(10):
        Ra_2[i]=gg[i]*math.cos(theta)*abs((T_i[i]-T_r[i]))*(L)**3
        h_r[i]=(k_air/L)*(0.68+((0.67*(Ra_2[i]**0.25))/(gggg)**(4/9)))

    #h_o(Ambient Side Natural Convection)%%%%%%%%%%%%%%%%%%%
    h_o=1.8+3.8*v_amb

    #Sky Temperature%%%%%
    T_sky=list(map(lambda x:x-6 ,T_amb))

    for i in range(10):
        alpha[i]=(t_g*((alpha_c*beta_c)+((alpha_tg)*(1-beta_c))))-(etha_c[i]*beta_c*alpha_c)
        
        #Radiative Heat Transfer Coefficient%%
        h_rad[i]=0.88*5.6704*(10**(-8))*(T_c[i]+T_sky[i])*(((T_sky[i])**2)+((T_c[i])**2))
        
        #Some More Heat Transfer Coefficients%%%
        u_b[i]=((L_tg/k_tg)+(1/h_air[i]))**-1
        u_T[i]=((L_g/k_g) +(1/(h_o)))**-1
        u_tair[i]=((1/u_b[i])+(1/u_T[i]))**-1
        u_bb[i]=((1/h_air[i])+(1/h_r[i])+(L_i/k_i))**-1
        u_L[i]=u_bb[i]+u_b[i]

    #Outlet Air Temperature%%%%%%%%%%
    for i in range(10):
        if forced==1:
            #Duct Mass Flow Rate is an input for Forced Convection%%%
            mdot[i]=p_2*A_2*v[i]        
        else:
            #Duct Mass Flow Rate Should be Calculated for Natural Convection%
            mdot[i]=((2*((A_2*p_2)**2)*A_pv*I[i]*etha_th*g*L*math.sin(theta))/(T_2[i]*C_p*(2+((f[i]*L)/D_h))))**(1/3)
        
        zzz[i]=((u_L[i]-((u_b[i]**2)/(u_T[i]+u_b[i]+h_rad[i])))*(b/(mdot[i]*C_p)))

        T_airout[i]=(((((u_bb[i])*T_r[i]+((u_b[i]*(u_T[i]+h_rad[i]))/(u_T[i]+u_b[i]+h_rad[i]))*T_amb[i]-6*((h_rad[i]*u_b[i])/(u_T[i]+u_b[i]+h_rad[i]))+((u_b[i])/(u_T[i]+u_b[i]+h_rad[i]))*I[i]*alpha[i]))/((u_L[i]-((u_b[i]**2)/(u_T[i]+u_b[i]+h_rad[i])))))*(1-math.exp(-zzz[i]*L)))+T_r[i]*math.exp(-zzz[i]*L)

    TT=list(map(lambda x:x-273 ,T_airout)) # C to K

    #Backsurface Temperature%%%%%
    for i in range(10):
        T_bs[i]=(((h_air[i]+((u_tg*u_b[i])/(u_T[i]+u_b[i]+h_rad[i])))*((T_airout[i]+T_r[i])/2)+((u_tg*u_T[i])/(u_T[i]+u_b[i]+h_rad[i]))*T_amb[i]+((u_tg*h_rad[i])/(u_T[i]+u_b[i]+h_rad[i]))*T_sky[i]+I[i]*alpha[i]*(u_tg/(u_T[i]+u_b[i]+h_rad[i])))/(u_tg+h_air[i]))

    TTT=list(map(lambda x:x-273 ,T_bs)) # C to K

    #Solar Cell Temperature%%%%%%%%%%%%%%%
    for i in range(10):
        T_c[i]=(u_b[i]*((T_airout[i]+T_r[i])/2)+u_T[i]*T_amb[i]+h_rad[i]*T_sky[i]+I[i]*alpha[i])/(u_T[i]+u_b[i]+h_rad[i])

    TTTT=list(map(lambda x:x-273 ,T_c)) # C to K

    #heat gained%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i in range(10):
        qdot[i]=mdot[i]*C_p*(T_airout[i]-T_r[i])

    #Total Efficiency%%%%%%%%%%%%%%%%%%%%
    for i in range(10):
        electrical11[i]=etha_c[i]*(1-(0.006*(T_c[i]-298)))
        thermal[i]=(((qdot[i]+electrical11[i]*b*L*I[i])/(I[i]*b*L)))*100

    #Insulation Temperature%%%%%
    for i in range(10):
        tinsulation[i]=T_airout[i]-(((u_bb[i]*L_i)/(k_i))*(T_airout[i]-T_r[i]))

    Tx=list(map(lambda x:x-273 ,tinsulation)) # C to K

    #Electrical Efficiency%%%%%%%%%%%%%%%%%%%%
    for i in range(10):
        electrical1[i]=etha_c[i]*(1-(0.006*(T_c[i]-298)))

    for i in range(10):
        T_bs=list(map(lambda x:x+273 ,TTT)) 
        T_c=list(map(lambda x:x+273 ,TTTT))
        T_2=list(map(lambda x:x+273 ,TT))
        T_i=list(map(lambda x:x+273 ,Tx))
        etha_c[i]=electrical1[i]
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T_amb[i]=T_amb[i]
        v[i]=v[i]
        T_r[i]=T_r[i]
        I[i]=I[i]
        T_sky[i]=T_sky[i]

#outputs:
#Please Change (TT) below to your desired variable
print(f"The Temperature changes during the day: {TT}")

