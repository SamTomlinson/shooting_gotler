%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 RK                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uses a fourth order runge kutta method to march backwards in 
% space from the far field to the near field giving a boundary value
% at the near field. Outputs nodes and eigenmodes

function [eta, v] = RK(a,b,deltaeta,bcs,gotler,baseT,baseTdash,...
    kshoot,eigval)

    % set up variables used in RK and preallocate for speed
    eta = a:deltaeta:b; n = length(eta);
    v = zeros(length(bcs),n); v(:,n) = bcs;
    % runge kutta from far boundary in towards zero
    for k = n-1:-1:1 
        % for base flow
        j=5*k;
        % intermediate steps
        k1 = gotler(eta(k+1),v(:,k+1),baseT(j-1),...
            baseTdash(j-1),kshoot,eigval);
        k2 = gotler(eta(k+1)-0.5*deltaeta,v(:,k+1)...
            -0.5*deltaeta.*k1',baseT(j-2),baseTdash(j-2),kshoot,eigval);
        k3 = gotler(eta(k+1)-0.5*deltaeta,v(:,k+1)...
            -0.5*deltaeta.*k2',baseT(j-3),baseTdash(j-3),kshoot,eigval);
        k4 = gotler(eta(k+1)-deltaeta,v(:,k+1)-...
            deltaeta.*k3',baseT(j-4),baseTdash(j-4),kshoot,eigval);
        % next step
        v(:,k)=v(:,k+1)-deltaeta*(k1'+2.*k2'+2.*k3'+k4')./6;    
    end
    
end
    
    
    