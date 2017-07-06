function [eta, v] = RK(a,b,deltaeta,bcs,gotler,baseT,baseTdash,baseU,kappa,...
    beta,khat,intbaseT)

    eta = a:deltaeta:b; 
    n = length(eta);
    v = zeros(length(bcs),n);
    v(:,1) = bcs;
    
    for k = 1:(n-1) 
        k1 = gotler(eta(k),v(:,k),baseT(k),baseTdash(k),baseU(k),kappa,beta,...
            khat,intbaseT);
        k2 = gotler(eta(k)+0.5*deltaeta,v(:,k)+0.5*deltaeta*k1',baseT(k),...
            baseTdash(k),baseU(k),kappa,beta,khat,intbaseT);
        k3 = gotler(eta(k)+0.5*deltaeta,v(:,k)+0.5*deltaeta*k2',baseT(k),...
            baseTdash(k),baseU(k),kappa,beta,khat,intbaseT);
        k4 = gotler(eta(k)+deltaeta,v(:,k)+deltaeta*k3',baseT(k),baseTdash(k),...
            baseU(k),kappa,beta,khat,intbaseT);
        v(:,k+1)=v(:,k)+deltaeta*(k1'+2*k2'+2*k3'+k4')/6;
    end