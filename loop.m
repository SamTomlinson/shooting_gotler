%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 loop                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Imporves accuracy of shoot, iterates with increasing accuracy 



%                                 Key                                 % 
%
% tol - tolerance of shooting accuracy
%
% v - v0 and v0dash array 
%
% gotler - function containing de for gotler
%
% deltaeta - step size
%
% tol - tolerance
%
% bcs - values of boundary conditions (2D vector)
%
% init - initial guess, if not specified given as [-5,10]
%
% a,b - two ends of the domain
%
% flow parameters - gamma (specific heat), Pr (prandtl), C (sutherlands
% constant), D (fitting parameter), etab (matching point for edge of
% adjustment region)
% 
% eigval - found eigenvalue for entered khat
%
% k - spanwise wavenumber
% 
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE



%                           Mostly plotting                           %


function [eigval,H1]=loop(eigval,khat,a,b,A,deltaeta,a1,gotler,baseT,...
    baseTdash,shoot1,tol)

% Initialise vectors
    
vec=[]; eigvec=[];
    
% Loop through different khat values 
    
for shoot1=eigval-tol:tol/10:eigval+tol
        
% Far field boudary condition 
        
    a1 = [exp(-khat*b), -khat*exp(-khat*b)];
    a2 = [exp(-khat*A^2/(3*a^3)), khat*A^2/(a^4)*exp(-khat*A^2/(3*a^3))];
        
    % Runge kutta inwards
        
    [~, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,khat,shoot1);
    %[~, F1] = RK(a,b,deltaeta,a2,gotler,baseT,baseTdash,khat,shoot1);
    
% Boundary condition constraints
        
    H1=F1(2,1)-((khat*A^2)/(a^4))*F1(1,1);
        
    H2=F1(2,end) + khat*F1(1,end);
   
% Vector of H error and ks
        
    vec=[H1,vec]; eigvec=[shoot1,eigvec];
    
end
   
% Calculate the crossing points
    
zerIdx=[];
for i=1:length(vec)-1
    if ((vec(i)>0 && vec(i+1)<0) || (vec(i)<0 && vec(i+1)>0))
    zerIdx(end+1)=i; % save index of zero-crossing
    end
end

% Improve accuracy 

eigs=eigvec(zerIdx);
vecs=vec(zerIdx);
eigvalright=eigvec(zerIdx(1)-1);
eigvalleft=eigvec(zerIdx(1)+1);
vecsright=vec(zerIdx(1)-1);
vecsleft=vec(zerIdx(1)+1);

% Iterate in 

while (vecsleft>1e-16)
eignew=(eigvalleft+eigvalright)/2;
vecnew=(vecsleft+vecsright)/2;
if vecnew<0
    eigvalright=eignew;
    vecsright=vecnew;
end
if vecnew>0
    eigvalleft=eignew;
    vecsleft=vecnew;
end
end

eigval=eigvalleft;
