function [t,H,Q,p] = classicalrungekutta4(fcn, tspan, q0,p0, n)
t0=tspan(1);
t_end=tspan(2);
h = (t_end-t0)/n;
t = t0:h:t_end;
H = zeros(1,length(t));
Q = zeros(1,length(t));
p = zeros(1,length(t));
H(1,1)=0;
Q(1,1) = q0; 
p(1,1) = p0;

A=[0,0,0,0;1/2,0,0,0;0,1/2,0,0;0,0,1,0];
c=[0,1/2,1/2,1]';
b=[1/6,1/3,1/3,1/6]';


for i = 2:length(t)
    k_1 = feval(fcn,t(i-1)+c(1)*h,Q(1,i-1));
    k_2 = feval(fcn,t(i-1)+c(2)*h,Q(1,i-1)+h*c(2)*p(1,i-1)/0.9953+(h^2)*(A(2,1)*k_1));
    k_3 = feval(fcn,t(i-1)+c(3)*h,Q(1,i-1)+h*c(3)*p(1,i-1)/0.9953+(h^2)*(A(3,2)*k_2));
    k_4 = feval(fcn,t(i-1)+c(4)*h,Q(1,i-1)+h*c(4)*p(1,i-1)/0.9953+(h^2)*(A(4,3)*k_3));

    Q(1,i) = Q(1,i-1)+ h*p(1,i-1)/0.9953+ (h^2)*(b(1)*k_1+b(2)*k_2+b(3)*k_3+b(4)*k_4);
    p(1,i) = p(1,i-1)+(h)*(b(1)*k_1+b(2)*k_2+b(3)*k_3+b(4)*k_4);

    H(1,i) = abs(p(1,i)^2/(2*0.9953) + 90.5*0.4814E-3*(1-exp(-1.814*(Q(1,i)-q0)))^2 - p0^2/(2*0.9953));
end
end
