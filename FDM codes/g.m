function [res] = g(y,mm,N)
res = zeros(mm,1);
res(1) = y(2) + ((2/(pi*y(3)))^(1/2))*y(1);
res(2) = y(mm*N+2) - ((2/(pi*y(mm*N+3)))^(1/2))*y(mm*N+1);
res(3) = y(4) + 2*((2/(pi*y(3)))^(1/2))*(y(3)-1);
res(4) = y(mm*N+4) - 2*((2/(pi*y(mm*N+3)))^(1/2))*(y(mm*N+3)-1);
end
