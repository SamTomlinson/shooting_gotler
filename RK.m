%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 RK                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Uses a fourth order runge kutta method to march backward/forward in 
% time to get solution at opposite boudary



%                                 Key                                 % 
%
% eta - grid points
%
% v - v0 and v0dash array 
%
% gotler - function containing de for gotler
%
% deltaeta - step size
%
% bcs - values of boundary conditions (2D vector)
%
% a,b - two ends of the domain
% 
% k - spanwise wavenumber
% 
% eigval - eigenvalue shoot

% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE



%                             Runge Kutta                               %

function [eta, v] = RK(a,b,deltaeta,bcs,gotler,baseT,baseTdash,...
    kshoot,eigval)

    % Set up variables used in RK preallocate for speed
    
    eta = a:deltaeta:b; 
    n = length(eta);
    v = zeros(length(bcs),n);
    v(:,n) = bcs;
    %v(:,1) = bcs;
    
    % RK from far boundary in towards zero
    
    for k = n-1:-1:1 
        j=5*k;
        k;
        size(baseT);
        k1 = gotler(eta(k+1),v(:,k+1),baseT(j-1),baseTdash(j-1),kshoot,eigval);

        k2 = gotler(eta(k+1)-0.5*deltaeta,v(:,k+1)-0.5*deltaeta.*k1',baseT(j-2),...
            baseTdash(j-2),kshoot,eigval);

        k3 = gotler(eta(k+1)-0.5*deltaeta,v(:,k+1)-0.5*deltaeta.*k2',baseT(j-3),...
            baseTdash(j-3),kshoot,eigval);

        k4 = gotler(eta(k+1)-deltaeta,v(:,k+1)-deltaeta.*k3',baseT(j-4),baseTdash(j-4),...
            kshoot,eigval);

        v(:,k)=v(:,k+1)-deltaeta*(k1'+2.*k2'+2.*k3'+k4')./6;
        
    end
    
    
    