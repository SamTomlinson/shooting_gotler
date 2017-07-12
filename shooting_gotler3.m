%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           shooting_gotler                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Shooting method for solving 1D gotler stability equation with 
% Dirichlet BCS. 4th order RK and bisection to ensure BCs




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
% beta - streamwise wavenumber
%
% k - spanwise wavenumber
% 
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE



%                         Output eigenfuntion                          %



function [eta, v,eigval] = shooting_gotler3(gotler,deltaeta,a,b,khat) 

    % Parameters and base flow should really be put into funtion 

    Pr=1; C=0.509; D=1; A=3*(1+C)/Pr;
    
    
    % Solve for the base flow 
    
    [~,baseT,baseTdash]=baseflow(C,Pr,D,1,deltaeta,a,b);

    tic; % Begin time
    
    % Initialise vectors
    
    vec=[]; eigvec=[];
    
    % Loop through different khat values 
    
    for shoot1=0.01:0.01:0.6
        
        % Far field boudary condition 
        
        a1 = [exp(-khat*b), -khat*exp(-khat*b)];
        
        % Near field condition
        
        %a2= [exp(-shoot1*A^2/(3*a^3)), (shoot1*A^2/(a^4))*exp(-shoot1*A^2/(3*a^3))];
        
        % Runge kutta inwards
        
        [~, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,khat,shoot1);     
    
        % Boundary condition constraints
        
        H1=F1(2,1)-((khat*A^2)/(a^4))*F1(1,1);
        
        H2=F1(2,end) + khat*F1(1,end);
   
        % Vector of H error and ks
        
        vec=[H1,vec];
        eigvec=[shoot1,eigvec];
    
    end
    
    % Plot H vs eig
    
%     figure('position', [0,0,800,800]); 
%     plot(eigvec,vec,'k-','LineWidth',2); 
%     set(gca,'Fontsize',20)
%     ylabel('Near field error, $H(\hat{k})$','Interpreter',...
%         'LaTex','Fontsize',40)
%     xlabel('Wavenumber, $\hat{k}$','Interpreter', 'LaTex','Fontsize',40)
%     xlim([0.1,0.5])
%     %ylim([-1,1])
%     grid on
%     hold off;
    
    % Calculate the corssing point
    
    
    zerIdx=[];
    for i=1:length(vec)-1
        if ((vec(i)>0 && vec(i+1)<0) || (vec(i)<0 && vec(i+1)>0))
        zerIdx(end+1)=i; % save index of zero-crossing
        end
    end
    
    eigs=eigvec(zerIdx);
    eigval=eigs(1);
    
    
    % Plotting of solutions 
    
    a1 = [exp(-khat*b), -khat*exp(-khat*b)];
    [eta, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,khat,eigval);     
    H1=F1(2,1)-((khat*A^2)/(a^4))*F1(1,1);
    H2=F1(2,length(F1(2,:))) + khat*F1(1,length(F1(2,:)));
    
    %if (abs(H1) > tol)
    %    error('Zero bc not satiafied')
    %end
    %if (abs(H2) > tol)
    %    error('Far field bc not satisfied')
    %end
    
    v=F1;
    
%     figure('position', [0,0,800,800]); 
%     plot(eta,v(1,:),'k-','LineWidth',2); hold on; 
%     plot(eta,v(2,:),'r-','LineWidth',2); 
%     set(gca,'Fontsize',20)
%     l1=legend('$v_0(\eta)$','$v_{0\eta}(\eta)$');
%     set(l1, 'Interpreter','LaTex','Fontsize',30);
%     ylabel('Vel. in the temp. adj. region $v_0$','Interpreter',...
%         'LaTex','Fontsize',40)
%     xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
%     xlim([a,b])
%     grid on
%     hold off;
    
%     figure('position', [0,0,800,800]); 
%     plot(eta,-baseTdash.*v(1,:)./baseT,'k-','LineWidth',2);
%     set(gca,'Fontsize',20)
%      ylabel('Temp. in the temp. adj. region $T_0$','Interpreter',...
%         'LaTex','Fontsize',40)
%     xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
%     xlim([a,b])
%     grid on
%     hold off;
%     toc
    
end