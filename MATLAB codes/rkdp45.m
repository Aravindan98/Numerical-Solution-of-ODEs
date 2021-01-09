function [t,H,q,p] = EmbeddedRK45(fcn, tspan, q0,p0, n)
% classical fourth order method
% Solve the initial value problem
% y’ = fcn(t,y), t0 <= t <= b, y(t0)=y0

M = 0.9953;
S = 1.814;
D = 90.5*0.4814E-3;

t0 = tspan(1);
t_end = tspan(2);
h = (t_end-t0)/n; % step size.
t = t0:h:t_end;
disp(t)
H = zeros(1,length(t));
q = zeros(1,length(t)); % solution vector initialization
p = zeros(1,length(t)); % solution vector initialization
H(1,1)=0;
q(1,1) = q0; % initial condition
p(1,1) = p0; % initial condition
% J. R. Dormand and P. J. Prince coeffcients
A = [0, 0, 0, 0, 0, 0, 0;
    1/5, 0, 0, 0, 0, 0, 0;
    3/40, 9/40, 0, 0, 0, 0, 0;
    44/45, -56/15, 32/9, 0, 0, 0, 0;
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0;
    9017/3168, -355/33, 46732/5247, 49/176 , -5103/18656, 0, 0;
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
    

b1 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]';
b2 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]';

c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1]';

for i = 2:length(t)
    k_1 = feval(fcn,t(i-1),q(1,i-1));
    k_2 = feval(fcn,t(i-1)+c(2)*h,q(1,i-1)+h*c(2)*p(1,i-1)/M+(h^2)*A(2,1)*k_1);
    k_3 = feval(fcn,t(i-1)+c(3)*h,q(1,i-1)+h*c(3)*p(1,i-1)/M+(h^2)*A(3,1)*k_1+(h^2)*A(3,2)*k_2);
    k_4 = feval(fcn,t(i-1)+c(4)*h,q(1,i-1)+h*c(4)*p(1,i-1)/M+(h^2)*A(4,1)*k_1+(h^2)*A(4,2)*k_2+(h^2)*A(4,3)*k_3);
    k_5 = feval(fcn,t(i-1)+c(5)*h,q(1,i-1)+h*c(5)*p(1,i-1)/M+(h^2)*A(5,1)*k_1+(h^2)*A(5,2)*k_2+(h^2)*A(5,3)*k_3+(h^2)*A(5,4)*k_4);
    k_6 = feval(fcn,t(i-1)+c(6)*h,q(1,i-1)+h*c(6)*p(1,i-1)/M+(h^2)*A(6,1)*k_1+(h^2)*A(6,2)*k_2+(h^2)*A(6,3)*k_3+(h^2)*A(6,4)*k_4+(h^2)*A(6,5)*k_5);
    k_7 = feval(fcn,t(i-1)+c(7)*h,q(1,i-1)+h*c(7)*p(1,i-1)/M+(h^2)*A(7,1)*k_1+(h^2)*A(7,2)*k_2+(h^2)*A(7,3)*k_3+(h^2)*A(7,4)*k_4+(h^2)*A(7,5)*k_5+(h^2)*A(7,6)*k_6);

    q(1,i) = q(1,i-1)+ h*p(1,i-1)/M+ (h^2)*b1(1)*k_1+(h^2)*b1(2)*k_2+(h^2)*b1(3)*k_3+(h^2)*b1(4)*k_4+(h^2)*b1(5)*k_5+(h^2)*b1(6)*k_6+(h^2)*b1(7)*k_7;
    p(1,i) = p(1,i-1)+h*b2(1)*k_1+h*b2(2)*k_2+h*b2(3)*k_3+h*b2(4)*k_4+h*b2(5)*k_5+h*b2(6)*k_6+h*b2(7)*k_7;
    H(1,i) = abs(p(1,i)*p(1,i)/(2*M) + D*(1-exp(-S*(q(1,i)-q0)))*(1-exp(-S*(q(1,i)-q0))) - p0*p0/(2*M));
 
end
end