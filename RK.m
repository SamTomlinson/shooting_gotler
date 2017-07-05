function [eta, y] = RK(a,b,h,con,gotler,baseT,baseTdash,baseU,kappa,...
    betag,sigma,intbaseT)
    eta = a:h:b; 
    n = length(eta);
    y = zeros(length(con),n);
    y(:,1) = con;
    
    for k = 1:(n-1) 
        k1 = gotler(eta(k),y(:,k),baseT(k),baseTdash(k),baseU(k),kappa,betag,...
            sigma,intbaseT);
        k2 = gotler(eta(k)+0.5*h,y(:,k)+0.5*h*k1',baseT(k),...
            baseTdash(k),baseU(k),kappa,betag,sigma,intbaseT);
        k3 = gotler(eta(k)+0.5*h,y(:,k)+0.5*h*k2',baseT(k),...
            baseTdash(k),baseU(k),kappa,betag,sigma,intbaseT);
        k4 = gotler(eta(k)+h,y(:,k)+h*k3',baseT(k),baseTdash(k),...
            baseU(k),kappa,betag,sigma,intbaseT);
        y(:,k+1)=y(:,k)+h*(k1'+2*k2'+2*k3'+k4')/6;
    end