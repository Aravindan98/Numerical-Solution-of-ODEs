N=[10,20,40,80,160,320,640,1280];
t_init=0;
y_init=1;
t_end=5.2;
lambda=20;
error=zeros(1,length(N));
for i=1:length(N)
    n=N(i);
    [t,y_n]=forward_euler(t_init,y_init,t_end,n,@fcn,lambda);
    y_exact=exp(-lambda*t);
    error(i)=max(abs(y_exact-y_n));
end
fig=figure('units','inches','position',[0,0,6,6*0.7]);
p1=loglog(N,error,'-ok','markersize',6,'linewidth',1);
xlim([10,1280])
xticks(N)
xlabel('Number of points','fontsize',13,'color','k','interpreter','latex')
ylabel('Error','fontsize',13,'color','k','interpreter','latex')
title('Forward Euler(convergence)','fontsize',13,'color','k','interpreter','latex');

