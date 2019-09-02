clear all;
close all;
clc;

c   = 2;
res = 0.05;

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
        if (phi<=0)
            M(i,j) = 1;
        end
        if (psi<=0)
            H(i,j) = 1;
        end
    end
end
Mat = H.*M;

L1  = (-pi/2)*(ones(1,100));
L2  = (pi/2)*(ones(1,100));
L3  = linspace(-pi,pi,100);

figure(1)
set(gca,'fontsize',20)
set(gcf,'color','white');
hold on
surf(thetaA,thetadotA,Mat')
shading interp
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$\dot{\theta}$','fontsize',30,'interpreter','latex');
zlabel('Level Sets');
xlim([-pi,pi])
ylim([-pi pi])
% colorbar
% plot3(L1,L3,0*L2,'-w','linewidth',2.5);
% plot3(L2,L3,0*L2,'-w','linewidth',2.5);
view(2)
res = 1;
thetaA = [-pi:res:pi];
thetadotA = [-pi:res:pi];
Tspan = [0:0.01:100];
dt = Tspan(2)-Tspan(1);
x = zeros(length(Tspan),2);
t = 0;
for i = 1 : length(thetaA)
    for j = 1 : length(thetadotA)
        x0 = [thetaA(i)  thetadotA(j)];
        phi      = ((thetaA(i)^2) - (pi/2)^2) + (2*k*thetaA(i)*thetadotA(j));
        psi      = (thetaA(i)^2) - (pi/2)^2;
        xdot = dynamics(t,x0);
%         quiver3(x0(1),x0(2),Mat(i,j),xdot(1),xdot(2),0,'color','white','linewidth',1.5);
        if(phi<0 && abs(x0(1))<pi/2)
%         if(phi>=0 || abs(x0(1))>pi/2)
            xdot = dynamics(t,x0);
            x(1,:) = x0;
            for m = 1  : length(Tspan)-1
                x(m+1,:) = x(m,:) + (dt*dynamics(Tspan(m),x(m,:)));
            end
            plot3(x(:,1),x(:,2),25*ones(1,length(Tspan)),'-r','linewidth',1.5)
            hold on
            plot3(x(1,1),x(1,2),25*ones(1,length(Tspan)),'*g','linewidth',1.5)
            drawnow
        end
    end
end



function xdot = dynamics(t,x)

global k
xdot     = [0 0];
theta    = x(1);
thetadot = x(2); 
phi      = ((theta^2) - (pi/2)^2) + (2*k*theta*thetadot);
u        = 0 ; 
eta      = 1 ; 

if(phi>=0)
    c1   = k*theta ; 
    c2   = (-2*theta*thetadot)-(k*thetadot^2)-eta ; 
    u    = c2/c1;
end

xdot(1) = thetadot ; 
xdot(2) = u ; 

end

