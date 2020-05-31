%%% MEEN 646 Module 4 %%%

clear all
clc
mdot_d=3.5;
T2=300;
p02=10^5;
pi_ref=3.0;
om=3000; %rad/s
eta_d=0.86;
Dm2=0.16; 
Rm2=0.08; 
D3=0.3;
R3=0.15;
hb2=0.1;
hb3=0.04;
alp2=pi/2; 
bet3=pi/2;
mdot=[2.5:0.01:4.5];
R=287;
k=1.3975; %subjected to change
cp=1009; %subjected to change

a0=-0.3844; a1=3.0222; a2=-1.7778;

U2=om*(Dm2/2); U3=om*(D3/2);
IR2=Dm2/2-hb2/2; OR2=Dm2/2+hb2/2;
A2=pi*(OR2^2-IR2^2); 
A3=pi*D3*hb3;

for i=1:length(mdot)
    %Efficiency parameters (given)
    eta(i)=a0+a1*(mdot(i)/mdot_d)+a2*(mdot(i)/mdot_d)^2;
    p_ratio(i)=eta(i)*pi_ref;
    %Inlet parameters
    syms a %a=V2
    a=solve(p02==mdot(i)*R*T2/(A2*a)+mdot(i)*a/(2*A2),a);
    b(i,1:2)=double(vpa(a));
    V2(i)=double(vpa(min(a)));
    rho2(i)=mdot(i)/(A2*V2(i));
    P2(i)=rho2(i)*R*T2;
    %Outlet parameters
    P3(i)=P2(i)*p_ratio(i);
    T3s=T2*(p_ratio(i)^((k-1)/k));
    T3(i)=T2+(T3s-T2)/eta(i); 
    rho3(i)=P3(i)/(R*T3(i));
    W3(i)=mdot(i)/(rho3(i)*A3);
    %Velocity calculations
    bet2(i)=atan(V2(i)/U2);
    bet2_deg(i)=bet2(i)*180/pi; %in degree
    W2(i)=sqrt(V2(i)^2+U2^2);
    %check(i)=acos(U3/V3(i));
    V3(i)=sqrt(W3(i)^2+U3^2);
    alp3(i)=pi/2+atan(U3/W3(i));
    alp3_deg(i)=alp3(i)*180/pi; %in degree
    %Additional calculations
    phi(i)=W3(i)/U3;
    ome(i)=om;
    %Power(i)=(abs(mdot(i)*cp*(T3(i)-T2)))/10^3;
    Power(i)=mdot(i)*(cp*(T3(i)-T2)+0.5*(V3(i)^2-V2(i)^2)); %unit in W
    Power1(i)=Power(i)/10^3; %unit in kW
    i=i+1;
end 

figure (1)
plot(mdot',p_ratio','r')
xlabel('Mass flow [kg/s]')
ylabel('Pressure Ratio \pi')
%xlim([2 4.5])
grid on

figure (2)
plot(mdot',Power1','r')
xlabel('Mass flow [kg/s]')
ylabel('Power [kW]')
%xlim([3 4])
grid on

figure (3)
plot(mdot',eta','r')
xlabel('Mass flow [kg/s]')
ylabel('\eta')
%xlim([2 4.5])
grid on

figure (4)
plot(mdot',phi','r')
xlabel('Mass flow [kg/s]')
ylabel('\phi')
%xlim([3 4])
grid on

figure (5)
plot(mdot',ome','r')
xlabel('Mass flow [kg/s]')
ylabel('Rotation speed \omega')
%xlim([3 4])
grid on

figure (6)
plot(mdot',bet2_deg','r')
xlabel('Mass flow [kg/s]')
ylabel('\beta_2')
%xlim([3 4])
grid on

figure (7)
plot(mdot',alp3_deg','r')
xlabel('Mass flow [kg/s]')
ylabel('\alpha_3')
%xlim([3 4])
grid on

figure (8)
plot(mdot',V2','r*'); hold on
plot(mdot',W2','g'); hold off
xlabel('Mass flow [kg/s]')
ylabel('V_2,W_2(m/s)')
legend('V_2','W_2')
%xlim([3 4])
grid on

figure (9)
plot(mdot',V3','r*'); hold on
plot(mdot',W3','g'); hold off
xlabel('Mass flow [kg/s]')
ylabel('V_3,W_3(m/s)')
legend('V_3','W_3')
%xlim([3 4])
grid on
