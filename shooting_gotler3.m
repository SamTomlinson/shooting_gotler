%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           shooting_gotler                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shooting method for solving 1D gotler stability equation with 
% Dirichlet asymptotic BCS. Uses a 4th order RK and bisection method
% to ensure BCs need very accurate shoot. Outputs nodes, eigenmodes
% and the shoot value

function [eta, v,eigval] = shooting_gotler3(gotler,deltaeta,a,b,khat) 

    % parameters and base flow should really be put into funtion 
    Pr=1; C=0.509; D=1; A=3*(1+C)/Pr;
    % solve for the base flow 
    [~,baseT,baseTdash,~,~]=baseflow(C,Pr,D,deltaeta,a,b);
    % initialise vectors
    vec=[]; eigvec=[];
    % loop through different khat values 
    for shoot1=0.01:0.001:0.6
        % boundary conditions 
        a1 = [exp(-khat*b), -khat*exp(-khat*b)];
        % runge kutta inwards
        [~, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,...
            khat,shoot1);
        % boundary condition constraint
        H1=F1(2,1)-((khat*A^2)/(a^4))*F1(1,1);
        % vector of H error and ks
        vec=[H1,vec]; eigvec=[shoot1,eigvec];
    end  

    % plot H vs eig  
    figure('position', [0,0,800,800]); 
    plot(eigvec,vec,'k-','LineWidth',2); 
    set(gca,'Fontsize',20)
    ylabel('Near field error, $H(\hat{k})$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('Eig. val., $\beta^2(G^*-Q)$','Interpreter', ...
        'LaTex','Fontsize',40)
    xlim([0.01,0.6])
    grid on
    hold off;  

    % calculate index of the crossing points  
    zerIdx=[];
    for i=1:length(vec)-1
        if ((vec(i)>0 && vec(i+1)<0) || (vec(i)<0 && vec(i+1)>0))
            zerIdx(end+1)=i;
        end
    end
    % improve accuracy 
    eigs=eigvec(zerIdx); eigval=eigs(1);
    diff=1; tol=0.01;
    while abs(diff>1e-16)
        eigvalold=eigval
        [eigval,~]=loop(eigval,khat,a,b,A,deltaeta,a1,gotler,baseT,...
        baseTdash,shoot1,tol);
        diff=abs(eigvalold-eigval);
        tol=tol/10;
    end      
    % calculation of eigenmodes with accuarate khat    
    a1 = [exp(-khat*b), -khat*exp(-khat*b)];
    [eta, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,khat,eigval);     
    v=F1;
    % find local maxima
    maxs=[];
    for j=2:length(v(1,:))-1
        if (v(1,j-1)<v(1,j) && v(1,j)>v(1,j+1))
            maxs=[j,maxs];
        end
    end
    % normalise
    v(1,:)=(1/v(1,maxs(1)))*v(1,:);
    % eliminate deviation
    [~,b1]=find(v(1,:) == 1);
    [~,b2]=min(v(1,2:b1));
    v(:,1:b2)=0;

    % plotting of eigenmodes (if running evvsk % out)
    figure('position', [0,0,800,800]); 
    plot(eta,v(1,:),'LineWidth',2); 
    set(gca,'Fontsize',20)
    ylabel('Vel. in the temp. adj. region $v_0$','Interpreter',...
         'LaTex','Fontsize',40)
    xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex',...
        'Fontsize',40)
    xlim([0.1,b])
    grid on
    hold off; 
    
end