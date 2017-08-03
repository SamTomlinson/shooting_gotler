%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 loop                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Imporves accuracy of shoot, iterates with increasing accuracy by 
% zooming into crossing region and decreasing the mesh spacing untill
% machine accuracy reached. Outputs the shoot value and eigenmodes.

function [eigval,H1]=loop(eigval,khat,a,b,A,deltaeta,a1,gotler,baseT,...
    baseTdash,shoot1,tol)

    % initialise vectors  
    vec=[]; eigvec=[];
    % loop through different khat values on refined grid  
    for shoot1=eigval-tol:tol/10:eigval+tol      
        % far field boudary condition    
        a1 = [exp(-khat*b), -khat*exp(-khat*b)];
        % runge kutta inwards
        [~, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,khat,shoot1);
        % boundary condition constraints    
        H1=F1(2,1)-((khat*A^2)/(a^4))*F1(1,1);
        % vector of H error and shoots
        vec=[H1,vec]; eigvec=[shoot1,eigvec];
    end
    % calculate the crossing points   
    zerIdx=[];
    for i=1:length(vec)-1
        if ((vec(i)>0 && vec(i+1)<0) || (vec(i)<0 && vec(i+1)>0))
            zerIdx(end+1)=i; 
        end
    end
    
    % improve accuracy by middling procedure
    eigs=eigvec(zerIdx);
    vecs=vec(zerIdx);
    eigval=eigval(1);
    eigvalright=eigvec(zerIdx(1)-1);
    eigvalleft=eigvec(zerIdx(1)+1);
    vecsright=vec(zerIdx(1)-1);
    vecsleft=vec(zerIdx(1)+1);
    % iterate in 
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
    
end
