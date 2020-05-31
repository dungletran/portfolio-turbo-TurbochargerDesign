%                     MEEN-646: Module 2                                     
%           Design of Subsonic Turbine Blade M2-T                                                        

% Objective:                                                                              
%    Develop a design software that enables you to generate subsonic turbine blades 
%                                                                                        
% Given Parameters                                                                        
%    - Generate a family of profile (alpha1=90, alpha2=160), (alpha1=45, alpha2=160)              
%    - Blade chord C                                                                                     
%                                                                                         
% Instruction:                                                                            
%    Input: alpha1, alpha2 & iZone 

clc;

iZone=3;

%input:

%stator blade
alpha1 = 120; %in degree
alpha2 = 40; %in degree
a=0;

%rotor blade
% alpha1 = 38.5; %in degree
% alpha2 = 165.1; %in degree
% a=1;

%convert to rad
alpha1 = (alpha1/180)*pi;
alpha2 = (alpha2/180)*pi;

Cax_ratio = 0.4;
C = 40; %chord
s = 50; %spacing
n_iter = 1000; %number of iteration for camberline

n_iter_b1 = 19;
n_iter_b2 = 665;
n_iter_b3 = 995; %for zone3 only
x_cam1=n_iter_b1/n_iter*C;
x_cam2=n_iter_b2/n_iter*C;
x_cam3=n_iter_b3/n_iter*C;%for zone3 only

%cascade stagger angle
gamma = atan(sin(alpha2)/(-Cax_ratio*sin(alpha1-alpha2)/sin(alpha1)+cos(alpha2))); 
C_ax = C*sin(gamma);


%define camber line equation:
x_p0 = 0; y_p0 = 0;
x_p2 = C; y_p2 = 0;

%determine P1 coordinates by consider triangle P0P1P3
%a_1 = 1/3*C/sin(alpha1)*sin(gamma); %P1P0 length
%b_1 = 1/3*C/sin(alpha1)*sin(pi-alpha1-gamma); %P1P3 length
%c_1 = 1/3*C; %P0P3 length
%p = (a_1+b_1+c_1)/2;
%area = sqrt(p*(p-a_1)*(p-b_1)*(p-c_1));

%y_p1 = 2*area/(1/3*C);
%x_p1 = y_p1/tan(pi-alpha1-gamma);

%determine P1 coordinates, formula given in the book (equation 10.40)
phi1=pi/2-alpha1+gamma;
phi2=pi/2+alpha2-gamma;
y_p1 = C*(cot(phi1)/(1+cot(phi1)/cot(phi2)));
x_p1 = C*(1/(1+cot(phi1)/cot(phi2)));



%Bezier Curve
for i=1:1:n_iter
    zeta(i) = i/n_iter;
    x_cam(i) = (1-zeta(i))^2*x_p0 + 2*(1-zeta(i))*zeta(i)*x_p1 + zeta(i)^2*x_p2;
    y_cam(i) = (1-zeta(i))^2*y_p0 + 2*(1-zeta(i))*zeta(i)*y_p1 + zeta(i)^2*y_p2;
    
    %camber line tangent angle
    v_cam(i)=atan((-2*(1-zeta(i))*y_p0 + 2*(1-2*zeta(i))*y_p1 + 2*zeta(i)*y_p2)/(-2*(1-zeta(i))*x_p0 + 2*(1-2*zeta(i))*x_p1 + 2*zeta(i)*x_p2));
    
    %for the 2nd camberline
    y_cam1(i) = y_cam(i)+s;
end

%blade thickness
if iZone == 1
    for i=1:1:n_iter
        x(i)=x_cam(i)/C;
        if (x_cam(i)<x_cam1)
            t(i)=C*(0.3419*x(i)^0.4929);%zone1
        elseif (x_cam(i)>x_cam1) && (x_cam(i)<x_cam2)
            t(i)=C*(-15.631*x(i)^6 + 38.563*x(i)^5 - 38.22*x(i)^4 + 19.934*x(i)^3 - 6.2802*x(i)^2 + 1.1333*x(i) + 0.0307);%zone1
        elseif (x_cam(i)>x_cam2)
            t(i)=C*(75.656*x(i)^6 - 375.15*x(i)^5 + 774.1*x(i)^4 - 850.22*x(i)^3 + 524.07*x(i)^2 - 172.08*x(i) + 23.628);%zone1
        end
        
        if t(i)<0
           t(i)=0;
        end
    end
elseif iZone == 2
    for i=1:1:n_iter
        x(i)=x_cam(i)/C;
        if (x_cam(i)<x_cam1)
            t(i)=C*(0.6128*x(i)^0.4937);%zone2
        elseif (x_cam(i)>x_cam1) && (x_cam(i)<x_cam2)
            t(i)=C*(-35.559*x(i)^6 + 83.97*x(i)^5 - 79.529*x(i)^4 + 39.519*x(i)^3 - 11.876*x(i)^2 + 2.0934*x(i) + 0.0531);%zone2
        elseif (x_cam(i)>x_cam2)
            t(i)=C*(93.702*x(i)^6 - 455.68*x(i)^5 + 921.82*x(i)^4 - 991.96*x(i)^3 + 598.52*x(i)^2 - 192.34*x(i) + 25.931);%zone2 
        end
        
        if t(i)<0
           t(i)=0;
        end
    end
elseif iZone == 3
    for i=1:1:n_iter
        x(i)=x_cam(i)/C;
        if (x_cam(i)<x_cam1)
            t(i)=C*(0.8232*x(i)^0.4941);%zone3
        elseif (x_cam(i)>x_cam1) && (x_cam(i)<x_cam2)
            t(i)=C*(-56.476*x(i)^6 + 129.12*x(i)^5 - 118.24*x(i)^4 + 56.666*x(i)^3 - 16.456*x(i)^2 + 2.8703*x(i) + 0.0696);%zone3
        elseif (x_cam(i)>x_cam2) && (x_cam(i)<x_cam3)
            t(i)=C*(65.209*x(i)^6 - 309.36*x(i)^5 + 610.82*x(i)^4 - 641.17*x(i)^3 + 376.88*x(i)^2 - 118.1*x(i) + 15.713);%zone3
        elseif (x_cam(i)>x_cam3)
            t(i)=C*(7.1679*x(i)^2 - 14.367*x(i) + 7.1992);%zone3
        end
        
        if t(i)<0
           t(i)=0;
        end
    end
end 
    


%suction side coordinate
for i=1:1:n_iter
    x_S(i) = x_cam(i) - (t(i)/3.5)*sin(v_cam(i));
    y_S(i) = y_cam(i) + (t(i)/3.5)*cos(v_cam(i));
    y_S1(i) = y_S(i) + s;
end

%pressure side coordinate
for i=1:1:n_iter
    x_P(i) = x_cam(i) + (t(i)/3.5)*sin(v_cam(i));
    %x_P_test(i) = x_cam(i) - (t(i)/2)*sin(v_cam(i));
    y_P(i) = y_cam(i) - (t(i)/3.5)*cos(v_cam(i));
    y_P1(i) = y_P(i) + s;
end

if a==0
ang=gamma+3*pi/2;
else 
    ang=gamma+pi/2;
end 
for i=1:length(x_S)
G(i,1:2)=[sin(ang), cos(ang);-cos(ang), sin(ang)]*[x_S(i);y_S1(i)]-[0, 0;-cos(ang), sin(ang)]*[x_S(1);y_S1(1)];
E(i,1:2)=[sin(ang), cos(ang);-cos(ang), sin(ang)]*[x_cam(i);y_cam1(i)]-[0, 0;-cos(ang), sin(ang)]*[x_cam(1);y_cam1(1)];
F(i,1:2)=[sin(ang), cos(ang);-cos(ang), sin(ang)]*[x_P(i);y_P1(i)]-[0, 0;-cos(ang), sin(ang)]*[x_P(1);y_P1(1)];

H(i,1:2)=[sin(ang), cos(ang);-cos(ang), sin(ang)]*[x_S(i);y_S(i)];
I(i,1:2)=[sin(ang), cos(ang);-cos(ang), sin(ang)]*[x_cam(i);y_cam(i)];
J(i,1:2)=[sin(ang), cos(ang);-cos(ang), sin(ang)]*[x_P(i);y_P(i)];
i=i+1;
end 

figure(2)
plot(G(:,1),G(:,2));
hold on
plot(F(:,1),F(:,2));
plot(E(:,1),E(:,2));
plot(H(:,1),H(:,2));
plot(I(:,1),I(:,2));
plot(J(:,1),J(:,2));
axis([-20 70 -55 35])
%plot the blade
figure(1);

plot(x_cam,y_cam,'g')
hold on
plot(x_S,y_S,'r')
hold on
plot(x_P,y_P,'b')
hold on



%plot(x_P_test,y_P,'y')
%hold on
%for i=1:1000:(n_iter)
%    th = 0:pi/50:2*pi;
%    xunit = t(i)/2 * cos(th) + x_cam(i);
%    yunit = t(i)/2 * sin(th) + y_cam(i);
%    h = plot(xunit, yunit,'r');
%    hold on
%end

plot(x_cam,y_cam1,'g')
hold on
plot(x_S,y_S1,'r')
hold on
plot(x_P,y_P1,'b')
hold on


gamma=gamma*180/pi;
x_S=x_S'; x_P=x_P';
y_P=y_P'; y_S=y_S';

%plot(x_P_test,y_P1,'y')
%hold on
%for i=1000:1000:(n_iter-1000)
%    th = 0:pi/50:2*pi;
%    xunit = t(i)/2 * cos(th) + x_cam(i);
%    yunit = t(i)/2 * sin(th) + y_cam1(i);
%    h = plot(xunit, yunit,'r');
%    hold on
%end

xlim([-4 44])
%ylim([-0.2 1])
ylim([-30 90])
axis([-10 45 -30 25])

