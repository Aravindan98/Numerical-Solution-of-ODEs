clf;
clc;
clear;
tspan=[0 500]; % set time interval

q0 = 1.4155;
M=0.9953;
p0 = (1.545/48.888)*M;
S = 1.814;
D = 90.5*0.4814E-3;

h=[2, 2.3684, 2.3685];
%[t,y,q,p]=crk4(@rk,tspan,q0,p0,1000);
%[t1,y1,q1,p1]=FE(@rk,tspan,q0,p0,1000,h(1));
[t2,y2,q2,p2]=symplectic(@rk,tspan,q0,p0,2000,2.2);
%[t3,y3,q3,p3]=leapfrog(@rk,tspan,q0,p0,2000,2.2);

fig1=figure('units','inches','position',[0,0,6,6*0.7]);
plt=plot(t2(1,:),y2(1,:),'-k');
%p=plot(t2(1,:),y1(1,:),'-r',t3(1,:),q3(1,:),'-b');
%p=plot(t1(1,:),-p1(1,:),'-k',t2(1,:),-p2(1,:),'-r',t3(1,:),-p3(1,:),'-b');
xlabel('$t$','fontsize',13,'color','k','Interpreter','latex');
ylabel('$H(t)-H(0)$','fontsize',13,'color','k','fontweight','bold','interpreter','latex');
title('Symplectic h=','fontsize',13,'color','k','interpreter','latex');
%legend([p(1) p(2) p(3)],{'Forward Euler','Symplectic', 'LeapFrog'},'interpreter','latex');
%legend([p(1) p(2)],{'Symplectic', 'LeapFrog'},'interpreter','latex');


