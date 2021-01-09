function dv_dt = f(~,q)
         dv_dt = -2*1.814*(90.5*0.4814E-3)*exp(-S*(q-1.41))*(1-exp(-S*(q-1.41)))/0.9953;
end
