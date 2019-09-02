clear all
close all
clc

theta = -1*pi:0.4:1*pi;
thetadot=-1*pi:0.4:1*pi;

t = 0 ;
figure(1)
set(gcf,'color','white');
xlabel('$\mathbf{\theta}$','fontsize',44,'interpreter','latex');
ylabel('$\dot{\theta}$','fontsize',44,'interpreter','latex');
grid on
hold on


for i = 1 : length(theta)
    for j = 1 : length(thetadot)
        x = [theta(i) ; thetadot(j)];
        u = 0 ;
        
        % For PD control use
        if (theta(i)>pi/2)
            u = -[1,1]*x ; 
        elseif(theta(i)<-pi/2)
            u = -[1,1]*x ; 
        end
        
        
        % For heuristic repulsion
        if (theta(i)>pi/2)
            u = -1 ; 
        elseif(theta(i)<-pi/2)
            u = +1 ; 
        end
        
        
        
        xdot = dynamics(t,x,u);
        quiver(theta(i),thetadot(j),xdot(1),xdot(2),'MaxHeadSize',0.5,'AutoScaleFactor',0.4,'color','blue','linewidth',2);
        hold on
    end
end
set(gca,'fontsize',34)      
xlim([-1.5*pi,1.5*pi]);
ylim([-1.5*pi,1.5*pi]);
axis square

L1  = (-pi/2)*(ones(1,100));
L2  = (pi/2)*(ones(1,100));
L3  = linspace(-2*pi,2*pi,100);
plot(L1,L3,'-k','linewidth',2.5);
plot(L2,L3,'-k','linewidth',2.5);



function xdot = dynamics(t,x,u)

theta = x(1);
thetadot = x(2); 

xdot(1) = thetadot ; 
xdot(2) = u ; 

end