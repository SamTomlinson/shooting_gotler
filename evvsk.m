% ev vs wavenumber

ev=[];
for k=0.1:0.1:8
[eta, v,eigval] = shooting_gotler3(@gotler,0.001,0.1,20,k);
ev=[eigval, ev];
k
end

plot(0.1:0.1:8,ev(end-2:-1:1))