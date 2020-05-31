%%% MEEN 646 Module 3 %%%

%Initialize constants
% Assumptions:
% (1) Constant diameter and constant cross section 
% (2) V1=V1ax, incoming velocity is completely axial
% (3) Perfect gas

clc;


mdot=3.5;  
P1=1.6*1000*101.235;
T1=800;
pratio=1.2; %Pressure ratio
phi=0.4; % Stage flow coefficient
r=0.5; % Degree of reaction
eta=0.85; %Isentropic efficiency
Dm=[0.2:0.05:0.45];
V3ax=[50:0.05:100]; 
k=1.354; % Ratio of specific heat
R=287; % Gas constant
cp1=1098;
cp3=1090;
%alpha1=(90/180)*pi;

%Calculated values
P3=P1/pratio; % Outlet pressure
T3s=T1/((pratio)^((k-1)/k)); %Isentropic temperature
T3=T1-eta*(T1-T3s); %Temperature at exit
rho1=P1/(R*T1); %Density at inlet
rho3=P3/(R*T3); %Density at outlet
h1=cp1*T1;
h3=cp3*T3;

j=1;


for i=1:length(Dm)
    for j=1:length(V3ax)
        V1ax(j)=V3ax(j)*rho3/rho1;
        hb(i,j)=mdot/(rho1*pi*Dm(i)*V1ax(j)); %in m
        hb1(i,j)=mdot/(rho1*pi*Dm(i)*V1ax(j))*1000; %in mm
        Vm3=V3ax(j);
        omega(i,j)=Vm3/(phi*Dm(i)/2); %in rad/s
        omega1(i,j)=omega(i,j)/(2*pi/60); %in rpm
        pw(i,j)=mdot*((h1-h3)+0.5*(V1ax(j)^2-V3ax(j)^2)); %in W
        pw1(i,j)=mdot*((h1-h3)+0.5*(V1ax(j)^2-V3ax(j)^2))/1000; %in kW
        %pw(i,j)=mdot*omega(i,j)*Dm(i)/2*(omega(i,j)*Dm(i)/2+2*tan(alpha1-pi/2)*Vm3)/1000;%in kW
        lambda(i,j)=pw(i,j)/(mdot*(omega(i,j)*Dm(i)/2)^2);
        %lambda(i,j)=1+2*tan(alpha1-pi/2)*Vm3/(omega(i,j)*Dm(i)/2);
          
        
        %flow angles (assume constant mean diameter and meridional velocity
        alpha2(i,j)=acot((1/phi)*(lambda(i,j)/2-r+1))*180/(pi); %in degree
        alpha3(i,j)=acot((1/phi)*(-lambda(i,j)/2-r+1))*180/(pi); %in degree
        beta2(i,j)=acot((1/phi)*(lambda(i,j)/2-r))*180/(pi); %in degree
        beta3(i,j)=acot((1/phi)*(-lambda(i,j)/2-r))*180/(pi); %in degree
        if alpha2(i,j) < 0
            alpha2(i,j)=alpha2(i,j)+180;
        end
        if alpha3(i,j) < 0
            alpha3(i,j)=alpha3(i,j)+180;
        end
        if beta2(i,j) < 0
            beta2(i,j)=beta2(i,j)+180;
        end
        if beta3(i,j) < 0
            beta3(i,j)=beta3(i,j)+180;
        end
        j=j+1;         
    end 
    i=i+1;
end 

figure (1)
plot(pw1(1,:),omega1(1,:),'r')
hold on
plot(pw1(2,:),omega1(2,:),'g')
hold on
plot(pw1(3,:),omega1(3,:),'b')
hold on
plot(pw1(4,:),omega1(4,:),'c')
hold on
plot(pw1(5,:),omega1(5,:),'y')
hold on
plot(pw1(6,:),omega1(6,:),'m')
hold on
xlabel('Power [kW]')
ylabel('Frequency [rpm]')
legend('Dm=.20','Dm=.25','Dm=.30','Dm=.35','Dm=.40','Dm=.45')
grid on

figure (2)
plot(pw1(1,:),hb1(1,:),'r')
hold on
plot(pw1(2,:),hb1(2,:),'g')
hold on
plot(pw1(3,:),hb1(3,:),'b')
hold on
plot(pw1(4,:),hb1(4,:),'c')
hold on
plot(pw1(5,:),hb1(5,:),'y')
hold on
plot(pw1(6,:),hb1(6,:),'m')
hold on
xlabel('Power [kW]')
ylabel('Blade Height [mm]')
legend('Dm=.20','Dm=.25','Dm=.30','Dm=.35','Dm=.40','Dm=.45')
grid on

figure (3)
plot(omega1(1,:),hb1(1,:),'r')
hold on
plot(omega1(2,:),hb1(2,:),'g')
hold on
plot(omega1(3,:),hb1(3,:),'b')
hold on
plot(omega1(4,:),hb1(4,:),'c')
hold on
plot(omega1(5,:),hb1(5,:),'y')
hold on
plot(omega1(6,:),hb1(6,:),'m')
hold on
xlabel('Frequency [rpm]')
ylabel('Blade Height [mm]')
legend('Dm=.20','Dm=.25','Dm=.30','Dm=.35','Dm=.40','Dm=.45')
grid on

figure (4)
plot(omega1(1,:),lambda(1,:),'r')
hold on
plot(omega1(2,:),lambda(2,:),'g')
hold on
plot(omega1(3,:),lambda(3,:),'b')
hold on
plot(omega1(4,:),lambda(4,:),'c')
hold on
plot(omega1(5,:),lambda(5,:),'y')
hold on
plot(omega1(6,:),lambda(6,:),'m')
hold on
xlabel('Frequency [rpm]')
ylabel('Specific Stage Load Coefficient \lambda')
legend('Dm=.20','Dm=.25','Dm=.30','Dm=.35','Dm=.40','Dm=.45')
grid on

figure (5)
plot(lambda(1,:),pw1(1,:),'r')
hold on
plot(lambda(2,:),pw1(2,:),'g')
hold on
plot(lambda(3,:),pw1(3,:),'b')
hold on
plot(lambda(4,:),pw1(4,:),'c')
hold on
plot(lambda(5,:),pw1(5,:),'y')
hold on
plot(lambda(6,:),pw1(6,:),'m')
hold on
xlabel('\lambda')
ylabel('Power [kW]')
legend('Dm=.20','Dm=.25','Dm=.30','Dm=.35','Dm=.40','Dm=.45')
grid on

figure (6)
plot(alpha2(1,:),lambda(1,:),'r')
hold on
plot(alpha3(2,:),lambda(2,:),'g')
hold on
plot(beta2(3,:),lambda(3,:),'b')
hold on
plot(beta3(4,:),lambda(4,:),'c')
hold on
ylabel('\lambda')
xlabel('\alpha_2,\beta_2,\alpha_3,\beta_3')
legend('\alpha_2','\alpha_3','\beta_2','\beta_3')
grid on
