function vecq = gotler(x,q,baseT,baseTdash,betag,Gstar,Q,sigma)
    vecq(1) = q(2);
    vecq(2) = q(2)*baseTdash./baseT ...
        + ((baseT^(-2)).*betag.*q(1) + ...
        (Gstar-Q)*(baseTdash*betag/sigma))*q(1);
    
    