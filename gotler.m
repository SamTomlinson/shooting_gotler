%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               gotler                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Encode gotler system




%                                 Key                                 % 
%
% eta - grid points
%
% q - derivatives v_0 and v_0'
%
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE
%
% khat - spanwise wavenumber
%
% ev - eigen value 
%



%                          Gotler system                               %



function vecq = gotler(~,q,baseT,baseTdash,khat,ev)
    
% Diff of q1 is q2 

vecq(1) = q(2);
    
% Diff of q2 is the rest of the system

vecq(2) =  (2.*q(2).*baseTdash)./(baseT) ...
        + (baseT.^2).*(khat.^2).*q(1) ...
        + (1/ev).*baseTdash.*(khat.^2).*q(1);
    
    