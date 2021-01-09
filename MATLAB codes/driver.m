clear;
clc;
close all;
t_initial=0;
y_initial=1;
t_end=3.5;
lambda=10;
[t1,y1]=forward_euler(t_initial,y_initial,t_end,100,@fcn,lambda);
[t2,y2]=forward_euler(t_initial,y_initial,t_end,20,@fcn,lambda);
y1_exact=exp(-lambda*t1);
y2_exact=exp(-lambda*t2);

fig1=figure('units','inches','position',[0,0,6,6*0.7]);
%p1=plot(t1,y1,'ok','markersize',6,'linewidth',1);
p1=plot(t1,y1,'ok',t2,y2,'*r','markersize',6,'linewidth',1);
hold on;
p2=plot(t1,y1_exact,'linewidth',2);
%set(p2,'Color',[1,0,0])
xlabel('$t$','fontsize',13,'color','k','Interpreter','latex');
ylabel('$\bf{Solution}$','fontsize',13,'color','k','fontweight','bold','interpreter','latex');
title('Forward Euler','fontsize',13,'color','k','interpreter','latex');
%legend([p1 p2],{'N=100', 'Exact Solution'},'interpreter','latex');
%legend([p1(1) p1(2) p2],{'N=100', 'N=10', 'Exact Solution'},'interpreter','latex');

%fig2=figure('units','inches','position',[0,0,6,6*0.7]);
%p1=semilogy(t1,abs(y1-y1_exact),'ok','markersize',6,'linewidth',1);
%set(p2,'Color',[1,0,0])
%xlabel('$t$','fontsize',13,'color','k','Interpreter','latex');
%ylabel('\bf{Error}','fontsize',13,'color','k','fontweight','bold','interpreter','latex');
%title('Forward Euler','fontsize',13,'color','k','interpreter','latex');
%legend(p1(1),'N=100','interpreter','latex');
