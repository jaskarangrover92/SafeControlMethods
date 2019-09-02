clear all
close all
clc

c = 2;
res = 0.01;

thetaA = [-pi:res:pi];
thetadotA = [-pi:res:pi];

k        = 2;

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
colorbar
view(2)
axis square
plot3(L1,L3,14*L2,'-k','linewidth',1.5);
plot3(L2,L3,14*L2,'-k','linewidth',1.5);

res = 0.2;
thetaA = [-pi:res:pi];
thetadotA = [-pi:res:pi];
Tspan = [0:0.05:100];
dt = Tspan(2)-Tspan(1);
x = zeros(length(Tspan),2);
t = 0;
for i = 1 : length(thetaA)
    for j = 1 : length(thetadotA)
        x0 = [thetaA(i)  thetadotA(j)];
        xdot = dynamics(t,x0);
        quiver(x0(1),x0(2),xdot(1),xdot(2),'MaxHeadSize',0.5,'AutoScaleFactor',0.4,'color','white','linewidth',2);

        x(1,:) = x0;
        for m = 1  : length(Tspan)-1
            x(m+1,:) = x(m,:) + (dt*dynamics(Tspan(m),x(m,:)));
        end
        x(end,:)
        plot3(x(end,1),x(end,2),25,'*m','linewidth',1.5)
        drawnow
    end
end



function xdot = dynamics(t,x)

xdot     = [0 0];
theta    = x(1);
thetadot = x(2); 
k        = 2 ; 
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

