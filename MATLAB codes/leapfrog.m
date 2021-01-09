function [t,error,q,p] = LF(fcn, tspan, q0,p0, n, h)
% classical fourth order method
% Solve the initial value problem
% y’ = fcn(t,y), t0 <= t <= b, y(t0)=y0

M = 0.9953;
S = 1.814;
D = 90.5*0.4814E-3;

t0 = tspan(1);
t_end = tspan(1)+(n-1)*h;
h = (t_end-t0)/n; % step size.
t = t0:h:t_end;
error = zeros(1,length(t));
q = zeros(1,length(t)); % solution vector initialization
q_half = zeros(1,length(t));
p = zeros(1,length(t)); % solution vector initialization
error(1,1)=0;
q(1,1) = q0; % initial condition
p(1,1) = p0; % initial condition
q_half(1,1) = q0 + (h/2)*(p0/M);


for i = 2:length(t)
    
	p(1,i) = p(1,i-1) + h*feval(fcn,0,q_half(1,i-1));
	q_half(1,i) = h*p(1,i)/M + q_half(1,i-1);
	q(1,i) = 0.5*(q_half(1,i)+q_half(1,i-1));
	
    error(1,i) = p(1,i)*p(1,i)/(2*M) + D*(1-exp(-S*(q(1,i)-q0)))*(1-exp(-S*(q(1,i)-q0))) - p0*p0/(2*M);
    error(1,i) = abs(error(1,i));
end
end