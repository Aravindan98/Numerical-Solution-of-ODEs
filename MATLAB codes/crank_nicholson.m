function[t,y]=crank_nicholson(t_initial,y_initial,t_end,n,fcn,lambda)
h=(t_end-t_initial)/n;
t=linspace(t_initial,t_end,n)';
y=zeros(n,1);
y(1)=y_initial;
for i=2:n
    y(i)=y(i-1)+h*feval(fcn,t(i-1),y(i-1),lambda);  %find forward estimate
    y(i)=y(i-1)+(h/2)*(feval(fcn,t(i-1),y(i-1),lambda)+feval(fcn,t(i),y(i),lambda))
end
end