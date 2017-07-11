function vecq = gotler(~,q,baseT,baseTdash,kshoot,eigval)
    vecq(1) = q(2);
    vecq(2) =  (2.*q(2).*baseTdash)./(baseT) ...
        + (baseT.^2).*(kshoot.^2).*q(1) ...
        + (1/eigval).*baseTdash.*(kshoot.^2).*q(1);
    
    