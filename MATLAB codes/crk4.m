function [t,error,q,p] = classicalrungekutta4(fcn, tspan, q0,p0, n)
t0=tspan(1);
t_end=tspan(2);
h = (t_end-t0)/n; % step size.
t = t0:h:t_end;
error = zeros(1,length(t));
q = zeros(1,length(t)); % solution vector initialization
p = zeros(1,length(t)); % solution vector initialization
error(1,1)=0;
q(1,1) = q0; % initial condition
p(1,1) = p0;

A=[0,0,0,0;1/2,0,0,0;0,1/2,0,0;0,0,1,0];
c=[0,1/2,1/2,1]';
b=[1/6,1/3,1/3,1/6]';

M = 0.9953;
S = 1.814;
D = 90.5*0.4814E-3;
for i = 2:length(t)
    k_1 = feval(fcn,t(i-1)+c(1)*h,q(1,i-1));
    k_2 = feval(fcn,t(i-1)+c(2)*h,q(1,i-1)+h*c(2)*p(1,i-1)/M+(h^2)*(A(2,1)*k_1));
    k_3 = feval(fcn,t(i-1)+c(3)*h,q(1,i-1)+h*c(3)*p(1,i-1)/M+(h^2)*(A(3,2)*k_2));
    k_4 = feval(fcn,t(i-1)+c(4)*h,q(1,i-1)+h*c(4)*p(1,i-1)/M+(h^2)*(A(4,3)*k_3));

    q(1,i) = q(1,i-1)+ h*p(1,i-1)/M+ (h^2)*(b(1)*k_1+b(2)*k_2+b(3)*k_3+b(4)*k_4);
    p(1,i) = p(1,i-1)+(h)*(b(1)*k_1+b(2)*k_2+b(3)*k_3+b(4)*k_4);
    
    error(1,i) = p(1,i)^2/(2*M) + D*(1-exp(-S*(q(1,i)-q0)))*(1-exp(-S*(q(1,i)-q0))) - p0*p0/(2*M);
    error(1,i) = abs(error(1,i));
end
end  