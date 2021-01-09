 D = 90.5*0.4814E-3;
 S = 1.814;
 q0 = 1.4155;
 M = 0.9953;
 v0 = 1.545 /48.888;
tspan = [0 100];
IC = [q0 v0];
opt = odeset('AbsTol',1E-4);
[t, Q] = ode45(@dqdtfn,tspan,IC,opt);
q = Q(:,1);
v = Q(:,2);
% Classically E would be (1/2)Mv^2 + U
% If that applies here, then:
E = 1/2*M*v.^2 + D*(1- exp(-S*(q-q0)).^2)-(1/2);
subplot(2,1,1)
plot(t,q),grid
xlabel('t'),ylabel('q')
subplot(2,1,2)
plot(t,E),grid
xlabel('t'),ylabel('E')
function dqdt = dqdtfn(~,Q)
         D = 90.5*0.4814E-3;
         S = 1.814;
         q0 = 1.41;
         M = 0.9953;
         
         q = Q(1);
         v = Q(2);
         dvdt = -2*D*S*exp(-S*(q-q0))*(1-exp(-S*(q-q0)))/M;
         dqdt = [v; dvdt];  
end