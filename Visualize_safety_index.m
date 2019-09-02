% This code plots a surface representation of the safety index in the state
% space of the 2D system considered in class phi = ((theta^2) - (pi/2)^2) + (2*k*theta*thetadot) ;
clear all;
close all;
clc;

c = 2;
res = 0.01;

thetaA = [-pi:res:pi];
thetadotA = [-pi:res:pi];

global k 
k = 1 ;

M   = zeros(length(thetaA),length(thetadotA));
H   = zeros(length(thetaA),length(thetadotA));

for  i = 1 : length(thetaA)
    for j = 1 : length(thetadotA)
        theta    = thetaA(i);
        thetadot = thetadotA(j) ;
        phi = ((theta^2) - (pi/2)^2) + (2*k*theta*thetadot) ;
        psi = (theta^2) - (pi/2)^2;
        M(i,j) = phi;
        
    end
end
Mat = M;

figure(1)
set(gca,'fontsize',20)
set(gcf,'color','white');
hold on
surf(thetaA,thetadotA,Mat')
shading interp
hold on
surf(thetaA,thetadotA,0*Mat')

xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$\dot{\theta}$','fontsize',30,'interpreter','latex');
zlabel('$\phi$','fontsize',30,'interpreter','latex');
title('$\mbox{Safety Index } \phi$','fontsize',30,'interpreter','latex');

xlim([-pi,pi])
ylim([-pi pi])
colorbar
grid on
view([45 45])
% view(2)
