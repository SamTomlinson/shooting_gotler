%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               gotler                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Encodes the gotler differential equation. Outputs the solution 
% vector of eigenmode at node and its derivative

function vecq = gotler(~,q,baseT,baseTdash,khat,ev)
    
    % differential of q1 is q2 
    vecq(1) = q(2);    
    % differential of q2 is the rest of the system
    vecq(2) =  (2.*q(2).*baseTdash)./(baseT) ...
        + (baseT.^2).*(khat.^2).*q(1) ...
        + (1/ev).*baseTdash.*(khat.^2).*q(1);
    
end
    
    