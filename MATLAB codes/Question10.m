tspan=[0 500]; 

q0 = 1.4155;
M=0.9953;
p0 = (1.545/48.888)*M;
S = 1.814;
D = 90.5*0.4814E-3;

N=[250,500,1000,2000,4000];
[t2,y2,q2,p2]=classicalrungekutta4(@f,tspan,q0,p0,N(5));
[t1,y1,q1,p1]=EmbeddedRK45(@f,tspan,q0,p0,N(5));

%fig1=figure('units','inches','position',[0,0,6,6*0.7]);
%subplot(1,2,1);
%plot(t1(1,:),y1(1,:),'-.b');
%xlabel('$t$','fontsize',10,'color','k','Interpreter','latex');
%ylabel('$q(t)$','fontsize',13,'color','k','fontweight','bold','interpreter','latex');
%ylabel('$|H(q(t),p(t))-H(q(0),p(0))|$','fontsize',10,'color','k','fontweight','bold','interpreter','latex');
%title('$RKDP45_{N=1000}$','fontsize',8,'color','k','interpreter','latex');



fig1=figure('units','inches','position',[0,0,6,6*0.7]);
subplot(1,2,1);
plot(t1(1,:),q1(1,:),'-.b');
xlabel('$t$','fontsize',10,'color','k','Interpreter','latex');
ylabel('$q(t)$','fontsize',13,'color','k','fontweight','bold','interpreter','latex');
%ylabel('$|H(q(t),p(t))-H(q(0),p(0))|$','fontsize',10,'color','k','fontweight','bold','interpreter','latex');
title('$RKDP45_{N=4000}$','fontsize',8,'color','k','interpreter','latex');

subplot(1,2,2);
plot(t2(1,:),q2(1,:),'-.r');
xlabel('$t$','fontsize',10,'color','k','Interpreter','latex');
ylabel('$p(t)$','fontsize',13,'color','k','fontweight','bold','interpreter','latex');
%ylabel('$|H(q(t),p(t))-H(q(0),p(0))|$','fontsize',10,'color','k','fontweight','bold','interpreter','latex');
title('$ClassicalRK4_{N=4000}$','fontsize',8,'color','k','interpreter','latex');