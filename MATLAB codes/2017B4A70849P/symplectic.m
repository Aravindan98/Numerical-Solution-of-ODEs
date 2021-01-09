function [t,H,q,p] = symplectic(fcn, tspan, q0, p0, n, h)
t0=tspan(1);
t_end = tspan(1)+(n-1)*h;
h = (t_end-t0)/n; 
t = t0:h:t_end;
q = zeros(1,length(t)); 
p = zeros(1,length(t));
H = zeros(1,length(t));
p(1,1) = p0;
H(1,1)=0;
q(1,1) = q0; 


for i = 2:length(t)
    p(1,i)=p(1,i-1)+h*feval(fcn,t(i-1),q(1,i-1));
    q(1,i)=q(1,i-1)+h*(p(1,i)/0.9953);
    H(1,i) = abs(p(1,i)^2/(2*0.9953) + 90.5*0.4814E-3*(1-exp(-1.814*(q(1,i)-q0)))^2 - p0^2/(2*0.9953));
end
end        