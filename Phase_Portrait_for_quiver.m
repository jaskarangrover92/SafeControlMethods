% This code plots instanteneous velocities of the system for different
% control methods i.e. safe set, barrier certificate, potential field and
% sliding mode control for safe control
clear all
close all
clc

theta    =   -1*pi:0.8:1*pi;
thetadot =  -1*pi:0.8:1*pi;
global k
k        = 2 ;
t        = 0 ;


figure(1)
set(gcf,'color','white');

for m = 1:4
    
    subplot(2,2,m)
    
    xlabel('$\mathbf{\theta}$','fontsize',34,'interpreter','latex');
    ylabel('$\dot{\theta}$','fontsize',34,'interpreter','latex');
    grid on
    hold on
    
    for i = 1 : length(theta)
        for j = 1 : length(thetadot)
            
            x = [theta(i) ; thetadot(j)];
            if (m==1)
                str = 'Potential Field Method';
                u = pfm(x);
            elseif (m==2)
                str = 'Barrier Certificate Method';
                u = bc(x);
            elseif (m==3)
                str = 'Safe Set';
                u = ss(x);
            elseif (m==4)
                str = 'Sliding Mode';
                u = smc(x);
            end
            
            xdot = dynamics(t,x,u);
            quiver(theta(i),thetadot(j),xdot(1),xdot(2),'MaxHeadSize',0.5,'AutoScaleFactor',0.4,'color','blue','linewidth',2);
            hold on
            axis equal
            xlim([-1.5*pi 1.5*pi]);
            ylim([-1.5*pi 1.5*pi]);
            drawnow
            
        end
    end
    set(gca,'fontsize',24)
    title(str,'fontsize',24,'interpreter','latex');
    hold off
    
    
end

function xdot = dynamics(t,x,u)

theta    = x(1);
thetadot = x(2); 
xdot(1)  = thetadot ; 
xdot(2)  = u ; 

end

function u = pfm(x)

global      k
theta    = x(1) ; 
thetadot = x(2);
phi      = ((theta^2) - (pi/2)^2) + (2*k*theta*thetadot) ;
psi      = (theta^2) - (pi/2)^2;
ucap     = -[.2,.2]*x;
u        = ucap;
if (psi>0)
    u = ucap - sin(theta);
end


end

function u = bc(x)

global      k
theta    = x(1) ; 
thetadot = x(2);
phi      = ((theta^2) - (pi/2)^2) + (2*k*theta*thetadot) ;
psi      = (theta^2) - (pi/2)^2;
ucap     = -[.2,.2]*x;
lambda   = 1;
%CHOP SMALL A VALUES
A         = 2*k*theta;
b         = (-2*theta*thetadot) -(k*thetadot^2) -(lambda*phi);
u         = quadprog(2*eye(1),-2*ucap,A,b);

tol = 1e-1;
if (norm(theta)<tol)
    u = ucap;
end

end


function u = ss(x)
global      k
theta    = x(1) ; 
thetadot = x(2);
phi      = ((theta^2) - (pi/2)^2) + (2*k*theta*thetadot) ;
psi      = (theta^2) - (pi/2)^2;
ucap     = -[.2,.2]*x;

eta        = 2 ; 
A          = 2*k*theta;

b          = (-2*theta*thetadot) -(k*thetadot^2) -(eta);
ustar      = quadprog(2*eye(1),-2*ucap,A,b);
u          = ucap ;

if (phi>=0)
    u    = ustar;
end



end

function u = smc(x)
global      k
theta    = x(1) ; 
thetadot = x(2);
phi      = ((theta^2) - (pi/2)^2) + (2*k*theta*thetadot) ;
psi      = (theta^2) - (pi/2)^2;
ucap     = -[.2,.2]*x;
u        = ucap;
if (phi>0)
    u    = ucap - 1*(2*k*theta);
end

end



