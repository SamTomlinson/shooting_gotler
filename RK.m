function [eta, v] = RK(a,b,deltaeta,bcs,gotler,baseT,baseTdash,...
    kshoot,eigval)

    eta = a:deltaeta:b; 
    n = length(eta);
    v = zeros(length(bcs),n);
    v(:,n) = bcs;
    %v(:,1) = bcs;
    
    for k = n-1:-1:1 
        
        size(eta);
        size(baseT);
        
        k1 = gotler(eta(k+1),v(:,k+1),baseT(k+1),baseTdash(k+1),kshoot,eigval);

        k2 = gotler(eta(k+1)-0.5*deltaeta,v(:,k+1)-0.5*deltaeta.*k1',baseT(k+1),...
            baseTdash(k+1),kshoot,eigval);

        k3 = gotler(eta(k+1)-0.5*deltaeta,v(:,k+1)-0.5*deltaeta.*k2',baseT(k+1),...
            baseTdash(k+1),kshoot,eigval);

        k4 = gotler(eta(k+1)-deltaeta,v(:,k+1)-deltaeta.*k3',baseT(k+1),baseTdash(k+1),...
            kshoot,eigval);

        v(:,k)=v(:,k+1)-deltaeta*(k1'+2.*k2'+2.*k3'+k4')./6;
        
    end

        
%     for k = 1:(n-1) 
%         k1 = gotler(eta(k),v(:,k),baseT(k),baseTdash(k),kshoot,eigval);
%         k2 = gotler(eta(k)+0.5*deltaeta,v(:,k)+0.5*deltaeta.*k1',baseT(k),...
%              baseTdash(k),kshoot,eigval);
%         k3 = gotler(eta(k)+0.5*deltaeta,v(:,k)+0.5*deltaeta.*k2',baseT(k),...
%              baseTdash(k),kshoot,eigval);
%         k4 = gotler(eta(k)+deltaeta,v(:,k)+deltaeta.*k3',baseT(k),baseTdash(k),...
%              kshoot,eigval);
%         v(:,k+1)=v(:,k)+deltaeta*(k1'+2*k2'+2*k3'+k4')/6;
%     end
    
 