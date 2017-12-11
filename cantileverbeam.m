clear all
close all

% figure print on/off
print = false;

%colormap definition
c = autumn(3);
set(groot,'defaultAxesColorOrder',c)


%dir = 'C:\Users\Rasmus\Dropbox\Apps\ShareLaTeX\Windmill Structural Mechanics\Lib';

%% Intialisation of Paramaters
r = [ 3.00 6.46 9.46 12.46 15.46 18.46 21.46 24.46 ...
    27.46 28.96 29.86 30.55 ];
l = r(2:end)-r(1:end-1);
EI1 = [ 1.7e9 5.6657e8 2.3568e8 1.1916e8 5.9832e7 2.9763e7 ...
    1.3795e7 5.3929e6 1.5133e6 6.3600e5 3.5733e5 1.5917e5 ];
EI2 = [ 1.7e9 1.4003e9 8.5164e8 5.2511e8 3.2937e8 2.0720e8 ...
    1.2099e8 5.9935e7 2.4543e7 1.4047e7 1.0060e7 7.2248e6 ];
m = [ 3.3e2 3.3719e2 2.7564e2 2.2902e2 1.9140e2 1.6692e2 ...
    1.5947e2 8.4519e1 4.7877e1 3.3029e1 2.5357e1 1.9990e1 ];
beta = [ 0.0 9.11 7.90 6.71 5.71 4.74 3.65 2.40...
    0.90 0.06 -0.44 -0.52 ];
pnorm = [2e3, 3e3, 1.2e3];
ptan = [170, 490, 200];
theta = [0, 0, 15.37];
V0 = [8,12,20];


%% Question 1
% Q1 is calculated for a beam with 
% the length of the blade 
% but with constant load and stiffness

% The deflection of a cantilever beam is calculated in the 
% matlab function cantBeamConstant, outputting the displace-
% ment field u(xi), and the tip angle alpha. The deflection
% is calculated using SoM#1 cantilever elementary case 3.



%Brug algoritmen fra opgave 2 til at beregne udb?jningen:

r1 = linspace(1,0);
pnorm1 = ones(1,length(r1));
ptan1 = zeros(1,length(r1));
EIy = ones(1,length(r1));
EIz = ones(1,length(r1));
theta1 = 0;
beta1 = zeros(1,length(r1));
[uz,uy] = cantBeamNumeric(pnorm1,ptan1,EIy,EIz,theta1,beta1,r1);

Q = 1;
EI = 1;
lb = 1;
[u,xi] = cantBeamAnalytic(Q,EI,lb);


figure(1)
clf
hold on

plot(xi,u,'linewidth',4);
plot(r1(10:10:end),uz(10:10:end),...
    'xk','linewidth',2, 'markersize',15);
legend('Given analytic solution','Coded solution')

grid on
%title('Cantilever beam bending')
xlabel('Position, \xi [-]')
ylabel('Beam deflection, [m]')
width=500;
height=200;
set(gcf,'units','points','position',[100,100,width,height])
box on
hold off
if print
    saveas(gca, fullfile(dir, 'cantileverbeam'), 'epsc');
end

unumeric = uz(10:10:end);
uanalytic = u(10:10:end);
error = uz(10:10:end)-u(10:10:end)./u(10:10:end);
numanalysis = [r1(10:10:end) ;...
    uz(10:10:end) ; u(10:10:end) ; error ];
%latex_table = latex(numanalysis);

%% Question 2
% Calculate the static deflection of the blade described 
% above using the static loads

[u_z,u_y] = cantBeamNumeric(pnorm,ptan,EI1,EI2,theta,beta,r);


figure(2)
clf
hold on

for i=1:length(pnorm)
    plot(r,u_z(i,:),'-x','markersize',10,'linewidth',2)
end
%title('Windmill-blade flapwise deflection')
legend('V_0 = 8','V_0 = 12','V_0 = 20', 'location','nw')
xlabel('Radial distance from nacelle, [m]')
ylabel('Deflection, [m]')
grid on
width=250;
height=350;
set(gcf,'units','points','position',[100,100,width,height])
box on
hold off

if print
    saveas(gca, fullfile(dir, 'normaldeflection'), 'epsc');
end

figure(3)
clf
hold on
for i=1:length(pnorm)
    plot(r,u_y(i,:),'-x','markersize',10,'linewidth',2)
end
%title('Windmill-blade tangential deflection')
legend('V_0 = 8','V_0 = 12','V_0 = 20', 'location','nw')
xlabel('Radial distance from nacelle, [m]')
ylabel('Deflection, [m]')
grid on
box on
width=250;
height=350;
set(gcf,'units','points','position',[100,100,width,height])

hold off
if print
    saveas(gca, fullfile(dir, 'tangentdeflection'), 'epsc');
end

%% Question 3
a = zeros(2*length(r));

for e = 1:2*length(r)
    
theta = 0;
    

p_1 = zeros(1,length(pnorm));
p_2 = zeros(1,length(pnorm));

u_y = zeros(length(theta),length(r));
u_z = zeros(length(theta),length(r));

for j=1:length(theta)
        p = zeros(1,2*length(r));
        p(e) = 1;
        
        for i = 1:length(r)
            pnorm(i) = p(2*i-1);
            ptan(i) = p(2*i);
        end
        
    for i=1:length(r)
        % Normal and tangential 
        % forces are transformed into the coordinates of
        % the shear-web spar, such that
        % the bending stiffnesses given can be
        % used.

        
        phi_p = beta(i)+theta(j);
        p_1(i) = (cosd(phi_p)*pnorm(i) + sind(phi_p)*ptan(i));
        p_2(i) = (cosd(phi_p)*ptan(i) - sind(phi_p)*pnorm(i));
    end
    
    T_y = zeros(1,length(r));
    T_z = zeros(1,length(r));
    M_y = zeros(1,length(r));
    M_z = zeros(1,length(r));
    M_1 = zeros(1,length(r));
    M_2 = zeros(1,length(r));
    kappa_1 = zeros(1,length(r));
    kappa_2 = zeros(1,length(r));
    kappa_y = zeros(1,length(r));
    kappa_z = zeros(1,length(r));
    delta_y = zeros(1,length(r));
    delta_z = zeros(1,length(r));
    
    
    for i=(length(r):-1:2)
        T_y(i-1) = T_y(i) + 1/2*(p_2(i-1)+p_2(i))*(r(i)-r(i-1));
        T_z(i-1) = T_z(i) + 1/2*(p_1(i-1)+p_1(i))*(r(i)-r(i-1));
        
        M_y(i-1) = M_y(i) - T_z(i)*(r(i)-r(i-1)) - ...
            (1/6*p_1(i-1)+1/3*p_1(i))*(r(i)-r(i-1))^2;
        M_z(i-1) = M_z(i) + T_y(i)*(r(i)-r(i-1)) + ...
            (1/6*p_2(i-1)+1/3*p_2(i))*(r(i)-r(i-1))^2;
    end
    
    for i=1:length(r)
        M_1(i) = M_y(i)*cosd(beta(i)+theta(j)) - ...
            M_z(i)*sind(beta(i)+theta(j));
        M_2(i) = M_y(i)*sind(beta(i)+theta(j)) + ...
            M_z(i)*cosd(beta(i)+theta(j));
        
        kappa_1(i) = M_1(i)/EI1(i);
        kappa_2(i) = M_2(i)/EI2(i);
        
        kappa_y(i) = kappa_1(i)*cosd(beta(i)+theta(j))...
            + kappa_2(i)*sind(beta(i)+theta(j));
        kappa_z(i) = -kappa_1(i)*sind(beta(i)+theta(j))...
            + kappa_2(i)*cosd(beta(i)+theta(j));
        
    end
    
    for i=(1:(length(r)-1))
        delta_y(i+1) = delta_y(i)+1/2*(kappa_y(i+1)...
            +kappa_y(i))*(r(i+1)-r(i));
        delta_z(i+1) = delta_z(i)+1/2*(kappa_z(i+1)...
            +kappa_z(i))*(r(i+1)-r(i));
        
        u_y(j,i+1) = u_y(j,i) + delta_z(i)*(r(i+1)-r(i))...
            + (1/6*kappa_z(1+i)+1/3*kappa_z(i))*(r(i+1)-r(i))^2;
        u_z(j,i+1) = u_z(j,i) - delta_y(i)*(r(i+1)-r(i))...
            - (1/6*kappa_y(1+i)+1/3*kappa_y(i))*(r(i+1)-r(i))^2;
    end
end
    
    for i=1:length(r)
        a(2*i-1,e) = u_z(i);
        a(2*i,e) = u_y(i);
    end
end

for i=1:length(r)
    M(i*2-1,i*2-1)=m(i);
    M(i*2,i*2)=m(i);
end

[V,D]=eigs(a*M,6);

omegaz1 = sqrt(1./D(1,1));
omegay1 = sqrt(1./D(2,2));
omegaz2 = sqrt(1./D(3,3));
omegay2 = sqrt(1./D(4,4));
omegaz3 = sqrt(1./D(5,5));
omegay3 = sqrt(1./D(6,6));


uz1 = zeros(length(r),1);
uz2 = uz1;
uz3 = uz1;
uy1 = uz1;
uy2 = uz1;
uy3 = uz1;
for i=1:length(r)
    uz1(i)=V(i*2-1,1);
    uz2(i)=V(i*2-1,3);
    uz3(i)=V(i*2-1,5);
    uy1(i)=V(i*2,2);
    uy2(i)=V(i*2,4);
    uy3(i)=V(i*2,6);
end


%Plotting Eigenmodes in the z-direction
figure(4)
clf
hold on

plot(r,uz1,'-x','linewidth',2,'markersize',10)
plot(r,uz2,'-x','linewidth',2,'markersize',10)
plot(r,uz3,'-x','linewidth',2,'markersize',10)

legend('EM(\omega_z1)','EM(\omega_z2)','EM(\omega_z3)',...
    'location','nw')
xlabel('Radial distance from nacelle, [m]')
ylabel('Deflection, [m]')
grid on
box on
width=500;
height=200;
set(gcf,'units','points','position',[100,100,width,height])

hold off


%Plotting Eigenmodes in the y-direction
figure(5)
clf
hold on

plot(r,uy1,'-x','linewidth',2,'markersize',10)
plot(r,uy2,'-x','linewidth',2,'markersize',10)
plot(r,uy3,'-x','linewidth',2,'markersize',10)

legend('EM(\omega_y1)','EM(\omega_y2)','EM(\omega_y3)',...
    'location','nw')
xlabel('Radial distance from nacelle, [m]')
ylabel('Deflection, [m]')
grid on
box on
width=500;
height=200;
set(gcf,'units','points','position',[100,100,width,height])

hold off


%% Question 4


% Flow specification:
theta = 1;
phi_0 = 8;
alpha_0 = phi_0-theta;
Vrel = 70;
t = linspace(0,pi);
x = sin(2*t);
x_dot = 2*cos(2*t);
y = sin(2*t);
y_dot = 2*cos(2*t);

% Airfoil data (Tjaere11.dat):
load tjaere11.dat
alphadata = tjaere11(:,1);
Cldata = tjaere11(:,2);
Cddata = tjaere11(:,3);

c = 1.5;

% Ambient conditions:
rho = 1.225;

phi_z = atand( (Vrel*sind(phi_0)-x_dot)./(Vrel*cosd(phi_0)) );
phi_y = atand( (Vrel*sind(phi_0))./(Vrel*cosd(phi_0)-y_dot) );

alpha_z = phi_z-theta;
alpha_y = phi_y-theta;


Cl_z = interp1(alphadata,Cldata,phi_z);
Cd_z = interp1(alphadata,Cddata,phi_z);

Cl_y = interp1(alphadata,Cldata,phi_y);
Cd_y = interp1(alphadata,Cldata,phi_y);

Cz = Cl_z.*cosd(phi_z)+Cd_z.*sind(phi_z);
Cy = Cl_y.*cosd(phi_y)+Cd_y.*sind(phi_y); 

Fz = 1/2*rho*Vrel^2*Cz;
Fy = 1/2*rho*Vrel^2*Cy;

Wz = trapz(t,Fz.*x_dot);
Wy = trapz(t,Fy.*y_dot);



%% Utilized functions:
function [u,xi] = cantBeamAnalytic(Q,EI,lb)
% Calculation of cantilever beam displacement distribution due to a
% concentrated load, from given analytic distribution.
x = linspace(0,1);
xi = (lb-x)./lb;
u = Q*lb^4/(24*EI)*(xi.^4-4*xi+3);
end

function [u_z,u_y] = cantBeamNumeric(pnorm,ptan,...
    EI1,EI2,theta,beta,r)

p_1 = zeros(1,length(pnorm));
p_2 = zeros(1,length(pnorm));

u_y = zeros(length(theta),length(r));
u_z = zeros(length(theta),length(r));
for j=1:length(theta)
    for i=(1:length(r))
        % Normal and tangential forces are 
        % transformed into the coordinates of
        % the shear-web spar, such that the bending 
        % stiffnesses given can be
        % used.
        phi_p = beta(i)+theta(j);
        p_1(i) = (cosd(phi_p)*pnorm(j) + sind(phi_p)*ptan(j));
        p_2(i) = (cosd(phi_p)*ptan(j) - sind(phi_p)*pnorm(j));
    end
    
    T_y = zeros(1,length(r));
    T_z = zeros(1,length(r));
    M_y = zeros(1,length(r));
    M_z = zeros(1,length(r));
    M_1 = zeros(1,length(r));
    M_2 = zeros(1,length(r));
    kappa_1 = zeros(1,length(r));
    kappa_2 = zeros(1,length(r));
    kappa_y = zeros(1,length(r));
    kappa_z = zeros(1,length(r));
    delta_y = zeros(1,length(r));
    delta_z = zeros(1,length(r));

    
    for i=(length(r):-1:2)
        T_y(i-1) = T_y(i) + 1/2*(p_2(i-1)+p_2(i))*(r(i)-r(i-1));
        T_z(i-1) = T_z(i) + 1/2*(p_1(i-1)+p_1(i))*(r(i)-r(i-1));
        
        M_y(i-1) = M_y(i) - T_z(i)*(r(i)-r(i-1)) - ...
            (1/6*p_1(i-1)+1/3*p_1(i))*(r(i)-r(i-1))^2;
        M_z(i-1) = M_z(i) + T_y(i)*(r(i)-r(i-1)) + ...
            (1/6*p_2(i-1)+1/3*p_2(i))*(r(i)-r(i-1))^2;
    end
    
    for i=1:length(r)
        M_1(i) = M_y(i)*cosd(beta(i)+theta(j)) - ...
            M_z(i)*sind(beta(i)+theta(j));
        M_2(i) = M_y(i)*sind(beta(i)+theta(j)) + ...
            M_z(i)*cosd(beta(i)+theta(j));
        
        kappa_1(i) = M_1(i)/EI1(i);
        kappa_2(i) = M_2(i)/EI2(i);
        
        kappa_y(i) = kappa_1(i)*cosd(beta(i)+theta(j)) + ...
            kappa_2(i)*sind(beta(i)+theta(j));
        kappa_z(i) = -kappa_1(i)*sind(beta(i)+theta(j)) + ...
            kappa_2(i)*cosd(beta(i)+theta(j));

    end
    
    for i=(1:(length(r)-1))
        delta_y(i+1) = delta_y(i)+1/2*...
            (kappa_y(i+1)+kappa_y(i))*(r(i+1)-r(i));
        delta_z(i+1) = delta_z(i)+1/2*...
            (kappa_z(i+1)+kappa_z(i))*(r(i+1)-r(i));
        
        u_y(j,i+1) = u_y(j,i) + delta_z(i)*(r(i+1)-r(i))...
            + (1/6*kappa_z(1+i)+1/3*kappa_z(i))*(r(i+1)-r(i))^2;
        u_z(j,i+1) = u_z(j,i) - delta_y(i)*(r(i+1)-r(i))...
            - (1/6*kappa_y(1+i)+1/3*kappa_y(i))*(r(i+1)-r(i))^2;
    end
end

end


function [u,xi,Theta] = cantBeamConstant(Q,EI,lb)
% Calculation of cantilever beam displacement distribution 
% due to a
% constant distributed load
xi = linspace(0,1);
u = Q*lb^4/(4*EI)*(xi.^2-2/3*xi.^3+1/6*xi.^4);
Theta = Q*lb^3/(6*EI)*180/pi;
end
function [u,xi] = cantBeamPoint(Q,EI,lb,alpha)
% Calculation of cantilever beam displacement distribution 
% due to a
% concentrated load
xi = linspace(0,1);
xi1 = xi(xi<=alpha);
xi2 = xi(xi>alpha);
u1 = Q*lb^3/(6*EI)*xi1.^2.*(3*alpha-xi1);
u2 = Q*lb^3/(6*EI)*alpha^2*(3*xi2-alpha);

u = [u1,u2];
end