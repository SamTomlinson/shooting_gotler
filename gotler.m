function vecq = gotler(~,q,baseT,baseTdash,baseU,kappa,~,~,~,sigma,...
    intbaseT)
    vecq(1) = q(2);
    vecq(2) = 2*q(2)*baseTdash./baseT ...
        + (0.5*kappa*(baseU^2)-intbaseT)*(baseTdash*(sigma^2)*q(1));
    
    