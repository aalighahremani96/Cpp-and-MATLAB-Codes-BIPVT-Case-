clear all
close all
clc
tic

%inputs%%%%%%%%%%%
loop=2;
forced=1;
nu=0.000015;
ro=1.200;
mu=1.51*10^-5;
v_amb=2;
g=9.8;
A_pv=1;
b=1;
C_p=1005;
p_2=1.184;
k_c=0.036;
k_g=1;
k_i=0.035;
k_tg=0.033;
k_air=0.025;
L=1;
L_c=300*10^(-6);
L_g=0.003;
L_i=0.05;
L_tg=0.0005;
alpha_c=0.8;
beta_c=0.9;
alpha_tg=0.5;
theta=0;
etha_c=0.2.*ones(1,10);
etha_th=0.5;
t_g=0.95;
pr=0.72;
D_h=9/100;
ductdepth=0.06;
A_2=ductdepth*b;
I=[300 400 550 650 658 600 560 420 250 110];
T_r=[32 33 35 36 42 44 42 41 41 39]+273;
T_amb=[29 30 32 35 38 40 41 41 40 39]+273;
v=[2.67 1.98 1.42 1.87 1.73 1.67 1.77 1.7 1.83 1.77];

%Backsurface, Cell Temperature, Outlet air Temperature Insulation Temperature Initial Guesses%%%%%%%
T_bs=[38 43 47 51 57 56 55 49 46 42]+273;
T_c=[38 41 48 50 55 55 53 51 47 42]+273;
T_2=273+[33 36 39 42 46 46 46 45 43 40];
T_i=273+30.*ones(1,10);

%Calculation%%%%%%%%%%%%
for w=1:loop
 
    
u_tg=(L_tg./k_tg).^-1;
   
%(Ra)/((H)^3(deltaT))%%
gg=107000000*ones(1,10,w);

%Ra, Gr, f%%%%%%%%%%%%%%%%%%
for uu= 1:10
    
    Ra(1,uu,w)=abs((T_bs(1,uu,w)-T_i(1,uu,w))).*gg(1,uu,w).*ductdepth.^3;
    Gr(1,uu,w)=Ra(1,uu,w)./pr;
    f(1,uu,w)=1.906*(Gr(1,uu,w)./pr).^(1/12);

%h_air(Duct Air)%%%%%%%%%%%%%%%%%%%%%%%
if forced==0
    
    Raa(1,uu,w)=abs((T_bs(1,uu,w)-T_i(1,uu,w))).*(gg(1,uu,w).*nu.*sin(theta)).*ro.*L.^3;
    
%Natural Convection%%%%%%
    h_air(1,uu,w)=(k_air./D_h).*0.0965.*Raa(1,uu,w).^0.29;

else
    
    Re(1,uu,w)=(p_2.*v(1,uu,w).*D_h)./(mu);
if Re(1,uu,w)>2300;
%Forced Convection(Turbulent)%%%%%%    
    h_air(1,uu,w)=(k_air./D_h).*0.023.*(Re(1,uu,w).^0.8).*(pr.^0.4);
else
%Forced Convection(Laminar)%%%%%%    
    h_air(1,uu,w)= (k_air./D_h).*2.98.*ones(1,10)
   
end
end

end



%h_r(Room Side Natural Convection)%%%%%%%%%%%%%%%%%%%
for uuu=1:10
    
    Ra_2(1,uuu,w)=gg(1,uuu,w).*cos(theta).*abs((T_i(1,uuu,w)-T_r(1,uuu,w))).*(L).^3;
    gggg=1+((0.492./pr).^(9/16));
    h_r(1,uuu,w)=(k_air./L)*(0.68+((0.67.*(Ra_2(1,uuu,w).^0.25))/(gggg).^(4/9)));

end


%h_o(Ambient Side Natural Convection)%%%%%%%%%%%%%%%%%%%
h_o=1.8+3.8.*v_amb;

%Sky Temperature%%%%%
T_sky=T_amb-6;


for u=1:10
    
    alpha(1,u,w)=(t_g.*((alpha_c.*beta_c)+((alpha_tg).*(1-beta_c))))-(etha_c(1,u,w).*beta_c.*alpha_c);
    
    %Radiative Heat Transfer Coefficient%%
    h_rad(1,u,w)=0.88*5.6704*(10^(-8)).*(T_c(1,u,w)+T_sky(1,u,w)).*(((T_sky(1,u,w)).^2)+((T_c(1,u,w)).^2));
    
    %Some More Heat Transfer Coefficients%%%
    
    u_b(1,u,w)=((L_tg./k_tg)+(1./h_air(1,u,w))).^-1;
    u_T(1,u,w)=((L_g./k_g) +(1./(h_o))).^-1;
    u_tair(1,u,w)=((1./u_b(1,u,w))+(1./u_T(1,u,w))).^-1;
    u_bb(1,u,w)=((1/h_air(1,u,w))+(1/h_r(1,u,w))+(L_i/k_i)).^-1;
    u_L(1,u,w)=u_bb(1,u,w)+u_b(1,u,w);
    
end






%Outlet Air Temperature%%%%%%%%%%


for t=1:1:10;  
    
    if forced==1
    %Duct Mass Flow Rate is an input for Forced Convection%%%
    mdot(1,t,w)=p_2.*A_2.*v(1,t,w);
    
    else
    %Duct Mass Flow Rate Should be Calculated for Natural Convection%
    mdot(1,t,w)=((2.*((A_2.*p_2).^2).*A_pv.*I(1,t,w).*etha_th.*g.*L.*sin(theta))./(T_2(1,t,w).*C_p.*(2+((f(1,t,w).*L)./D_h)))).^(1/3);
    
    end
    
    zzz(1,t,w)=((u_L(1,t,w)-((u_b(1,t,w).^2)./(u_T(1,t,w)+u_b(1,t,w)+h_rad(1,t,w)))).*(b./(mdot(1,t,w).*C_p)));

    T_airout(1,t,w)=(((((u_bb(1,t,w)).*T_r(1,t,w)+((u_b(1,t,w).*(u_T(1,t,w)+h_rad(1,t,w)))./(u_T(1,t,w)+u_b(1,t,w)+h_rad(1,t,w))).*T_amb(1,t,w)-6*((h_rad(1,t,w).*u_b(1,t,w))./(u_T(1,t,w)+u_b(1,t,w)+h_rad(1,t,w)))+((u_b(1,t,w))./(u_T(1,t,w)+u_b(1,t,w)+h_rad(1,t,w))).*I(1,t,w).*alpha(1,t,w)))./((u_L(1,t,w)-((u_b(1,t,w).^2)./(u_T(1,t,w)+u_b(1,t,w)+h_rad(1,t,w)))))).*(1-exp(-zzz(1,t,w).*L)))+T_r(1,t,w).*exp(-zzz(1,t,w).*L);
end
TT=T_airout-273;


%Backsurface Temperature%%%%%
for ttt=1:1:10;
    T_bs(1,ttt,w)=(((h_air(1,ttt,w)+((u_tg.*u_b(1,ttt,w))./(u_T(1,ttt,w)+u_b(1,ttt,w)+h_rad(1,ttt,w)))).*((T_airout(1,ttt,w)+T_r(1,ttt,w))./2)+((u_tg.*u_T(1,ttt,w))./(u_T(1,ttt,w)+u_b(1,ttt,w)+h_rad(1,ttt,w))).*T_amb(1,ttt,w)+((u_tg.*h_rad(1,ttt,w))./(u_T(1,ttt,w)+u_b(1,ttt,w)+h_rad(1,ttt,w))).*T_sky(1,ttt,w)+I(1,ttt,w).*alpha(1,ttt,w).*(u_tg./(u_T(1,ttt,w)+u_b(1,ttt,w)+h_rad(1,ttt,w))))./(u_tg+h_air(1,ttt,w)));
end

TTT=T_bs-273;

%Solar Cell Temperature%%%%%%%%%%%%%%%
for tttt=1:1:10;
    T_c(1,tttt,w)=(u_b(1,tttt,w).*((T_airout(1,tttt,w)+T_r(1,tttt,w))./2)+u_T(1,tttt,w).*T_amb(1,tttt,w)+h_rad(1,tttt,w).*T_sky(1,tttt,w)+I(1,tttt,w).*alpha(1,tttt,w))./(u_T(1,tttt,w)+u_b(1,tttt,w)+h_rad(1,tttt,w));
end

TTTT=T_c-273;

%heat gained%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ttttt=1:10;
    qdot(1,ttttt,w)=mdot(1,ttttt,w).*C_p.*(T_airout(1,ttttt,w)-T_r(1,ttttt,w));
end

%Total Efficiency%%%%%%%%%%%%%%%%%%%%
for tttttt=1:10
    electrical11(1,tttttt,w)=etha_c(1,tttttt,w).*(1-(0.006.*(T_c(1,tttttt,w)-298)));
    thermal(1,tttttt,w)=(((qdot(1,tttttt,w)+electrical11(1,tttttt,w).*b.*L.*I(1,tttttt,w))./(I(1,tttttt,w).*b.*L)))*100;
end


%Insulation Temperature%%%%%
for ttttttt=1:10
    tinsulation(1,ttttttt,w)=T_airout(1,ttttttt,w)-(((u_bb(1,ttttttt,w).*L_i)./(k_i)).*(T_airout(1,ttttttt,w)-T_r(1,ttttttt,w)));
end
TX=tinsulation-273;


%Electrical Efficiency%%%%%%%%%%%%%%%%%%%%
for tttttttt=1:10
    electrical1(1,tttttttt,w)=etha_c(1,tttttttt,w).*(1-(0.006.*(T_c(1,tttttttt,w)-298)));
end


for mm=1:10
    
    T_bs(1,mm,w+1)=TTT(1,mm,w)+273;
    T_c(1,mm,w+1)=TTTT(1,mm,w)+273;
    T_2(1,mm,w+1)=TT(1,mm,w)+273;
    T_i(1,mm,w+1)=TX(1,mm,w)+273;
    etha_c(1,mm,w+1)=electrical1(1,mm,w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_amb(1,mm,w+1)=T_amb(1,mm,w);
    v(1,mm,w+1)=v(1,mm,w);
    T_r(1,mm,w+1)=T_r(1,mm,w);
	I(1,mm,w+1)=I(1,mm,w);
    T_sky(1,mm,w+1)=T_sky(1,mm,w);

end
end

%Figures%%%%%%%%%%%%
x=1:10;
figure(1)
plot(x,TTTT(:,:,w),'*-k',x,TTT(:,:,w),'v-b',x,TT(:,:,w),'x-r')
grid on
set(gca,'fontname','Calibri (Body)','fontsize',12)
xticks([1:1:10])
yticks([0:5:65])
axis([1,10,0,65])
xlabel('Time(hour)')
ylabel('Temperature(°C)')
legend('Cell Temperature' ,'back surface Temperature','Outlet Air Temperature Temperature', 'Location', 'SouthEast')
xticklabels({'8AM','9AM','10AM','11AM','12PM','1PM','2PM','3PM','4PM','5PM'})

figure(2)
yyaxis left;
plot(x, electrical1(:,:,w)*100, '-o', 'LineWidth', 1.5);
ylabel('Electrical Efficiency');
xticklabels({'8AM','9AM','10AM','11AM','12PM','1PM','2PM','3PM','4PM','5PM'})
yticks([0:10:100])
axis([1,10,0,100])
yyaxis right;
plot(x, thermal(:,:,w), '--x', 'LineWidth', 1.5);
ylabel('Total Efficiency');
set(gca,'fontname','Calibri (Body)','fontsize',12)
grid on
yticks([0:10:100])
axis([1,10,0,100])
legend('Electrical Efficiency', 'Thermal Efficiency', 'Location', 'best');
toc

