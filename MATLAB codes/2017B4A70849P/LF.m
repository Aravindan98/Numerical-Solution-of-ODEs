function [t,H,q,p] = LF(fcn, tspan, q0,p0, n, h)

t0 = tspan(1);
t_end = tspan(1)+(n-1)*h;
h = (t_end-t0)/n;
t = t0:h:t_end;
H = zeros(1,length(t));
q = zeros(1,length(t));
q_half = zeros(1,length(t));
p = zeros(1,length(t));
H(1,1)=0;
q(1,1) = q0;
p(1,1) = p0; 
q_half(1,1) = q0 + (h/2)*(p0/0.9953);


for i = 2:length(t)

	p(1,i) = p(1,i-1) + h*feval(fcn,0,q_half(1,i-1));
	q_half(1,i) = h*p(1,i)/0.9953 + q_half(1,i-1);
	q(1,i) = 0.5*(q_half(1,i)+q_half(1,i-1));

    H(1,i) = abs(p(1,i)^2/(2*0.9953) + 90.5*0.4814E-3*(1-exp(-1.814*(q(1,i)-q0)))^2 - p0^2/(2*0.9953));
end
end
