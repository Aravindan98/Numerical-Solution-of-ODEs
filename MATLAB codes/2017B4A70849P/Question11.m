clc;
clear;
clf;
tspan=[0 500]; % set time interval

q_initial = 1.4155;
M=0.9953;
p_initial = (1.545/48.888)*M;
S = 1.814;
D = 90.5*0.4814E-3;

h=[2, 2.3684,2.3685, 2.2];

[t21,y21,q21,p21]=symplectic(@f,tspan,q_initial,p_initial,1000,h(1));
[t22,y22,q22,p22]=symplectic(@f,tspan,q_initial,p_initial,1000,h(2));
[t23,y23,q23,p23]=symplectic(@f,tspan,q_initial,p_initial,1000,h(3));
[t24,y24,q24,p24]=symplectic(@f,tspan,q_initial,p_initial,1000,h(4));

[t11,y11,q11,p11]=FE(@f,tspan,q_initial,p_initial,1000,h(1));
[t12,y12,q12,p12]=FE(@f,tspan,q_initial,p_initial,1000,h(2));
[t13,y13,q13,p13]=FE(@f,tspan,q_initial,p_initial,1000,h(3));
[t14,y14,q14,p14]=FE(@f,tspan,q_initial,p_initial,1000,h(4));

fig1=figure('units','inches','position',[0,0,6,6*0.7]);
subplot(2,2,1);
plot(t11(1,:),y11(1,:),'-b');
xlabel('$t$','fontsize',5,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',5,'color','k','fontweight','bold','interpreter','latex');
title('Forward Euler discrepancy h=2','fontsize',13,'color','k','interpreter','latex');
subplot(2,2,2);
plot(t12(1,:),y12(1,:),'-b');
xlabel('$t$','fontsize',5,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',5,'color','k','fontweight','bold','interpreter','latex');
title('Forward Euler discrepancy h=2.3684','fontsize',13,'color','k','interpreter','latex');
subplot(2,2,3);
plot(t13(1,:),y13(1,:),'-b');
xlabel('$t$','fontsize',13,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',7,'color','k','fontweight','bold','interpreter','latex');
title('Forward Euler discrepancy h=2.3685','fontsize',13,'color','k','interpreter','latex');
subplot(2,2,4);
plot(t14(1,:),y14(1,:),'-b');
xlabel('$t$','fontsize',13,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',7,'color','k','fontweight','bold','interpreter','latex');
title('Forward Euler discrepancy h=2.2','fontsize',13,'color','k','interpreter','latex');

fig2=figure('units','inches','position',[0,0,6,6*0.7]);
subplot(2,2,1);
plot(t21(1,:),y21(1,:),'-b');
xlabel('$t$','fontsize',5,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',5,'color','k','fontweight','bold','interpreter','latex');
title('Symplectic Euler discrepancy h=2','fontsize',13,'color','k','interpreter','latex');
subplot(2,2,2);
plot(t22(1,:),y22(1,:),'-b');
xlabel('$t$','fontsize',5,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',5,'color','k','fontweight','bold','interpreter','latex');
title('Symplectic Euler discrepancy h=2.3684','fontsize',13,'color','k','interpreter','latex');
subplot(2,2,3);
plot(t23(1,:),y23(1,:),'-b');
xlabel('$t$','fontsize',13,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',7,'color','k','fontweight','bold','interpreter','latex');
title('Symplectic Euler  discrepancy h=2.3685','fontsize',13,'color','k','interpreter','latex');
subplot(2,2,4);
plot(t24(1,:),y24(1,:),'-b');
xlabel('$t$','fontsize',13,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',7,'color','k','fontweight','bold','interpreter','latex');
title('Symplectic Euler  discrepancy h=2.2','fontsize',13,'color','k','interpreter','latex');

[t31,y31,q31,p31]=LF(@f,tspan,q_initial,p_initial,1000,h(1));
[t32,y32,q32,p32]=LF(@f,tspan,q_initial,p_initial,1000,h(2));
[t33,y33,q33,p33]=LF(@f,tspan,q_initial,p_initial,1000,h(3));
[t34,y34,q34,p34]=LF(@f,tspan,q_initial,p_initial,1000,h(4));

fig3=figure('units','inches','position',[0,0,6,6*0.7]);
subplot(2,2,1);
plot(t31(1,:),y31(1,:),'-b');
xlabel('$t$','fontsize',5,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',5,'color','k','fontweight','bold','interpreter','latex');
title('Leapfrog method discrepancy h=2','fontsize',13,'color','k','interpreter','latex');
subplot(2,2,2);
plot(t32(1,:),y32(1,:),'-b');
xlabel('$t$','fontsize',5,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',5,'color','k','fontweight','bold','interpreter','latex');
title('Leapfrog method discrepancy h=2.3684','fontsize',13,'color','k','interpreter','latex');
subplot(2,2,3);
plot(t33(1,:),y33(1,:),'-b');
xlabel('$t$','fontsize',13,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',7,'color','k','fontweight','bold','interpreter','latex');
title('Leapfrog method  discrepancy h=2.3685','fontsize',13,'color','k','interpreter','latex');
subplot(2,2,4);
plot(t34(1,:),y34(1,:),'-b');
xlabel('$t$','fontsize',13,'color','k','Interpreter','latex');
ylabel('$H(q(t),p(t))-H(q(0),p(0))$','fontsize',7,'color','k','fontweight','bold','interpreter','latex');
title('Leapfrog method discrepancy h=2.2','fontsize',13,'color','k','interpreter','latex');
