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



function [eta, v,kval] = shooting_gotler2(gotler,deltaeta,tol,a,b,eigval,init) 

    % Parameters and base flow should really be put into funtion 

    Pr=1; C=0.509; D=1; A=3*(1+C)/Pr;
    
    
    % Solve for the base flow 
    
    [~,baseT,baseTdash]=baseflow(C,Pr,D,1,deltaeta,a,b);

    tic; % Begin time
    
    % Initialise vectors
    
    vec=[]; kvec=[];
    
    % Loop through different khat values 
    
    for shoot1=0.1:0.001:8
        
        shoot1;
        % Far field boudary condition 
        
        a1 = [exp(-shoot1*b), -shoot1*exp(-shoot1*b)];
        
        % Runge kutta
        
        [~, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,shoot1,eigval);     
    
        % Boundary condition constraints
        
        H1=F1(2,1)-((shoot1*A^2)/(a^4))*F1(1,1);
        H2=F1(2,end) + shoot1*F1(1,end);
   
        % Vector of H error and ks
        
        vec=[H1,vec];
        kvec=[shoot1,kvec];
    
    end
    
    % Plot H vs k
    
    figure('position', [0,0,800,800]); 
    plot(kvec,vec,'k-','LineWidth',2); 
    set(gca,'Fontsize',20)
    ylabel('Near field error, $H(\hat{k})$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('Wavenumber, $\hat{k}$','Interpreter', 'LaTex','Fontsize',40)
    xlim([0.1,8])
    grid on
    hold off;
    
    % Calculate the corssing point
    
    dif = abs(vec-0);
    match = dif == min(dif);
    idx = find(dif == min(dif));
    closest = vec(idx);
    kval= kvec(idx);
    
    % Plotting of solutions 
    
    a1 = [exp(-kval*b), -kval*exp(-kval*b)];
    [eta, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,kval,eigval);     
    H1=F1(2,1)-((kval*A^2)/(a^4))*F1(1,1);
    H2=F1(2,length(F1(2,:))) + kval*F1(1,length(F1(2,:)));
    
    %if (abs(H1) > tol)
    %    error('Zero bc not satiafied')
    %end
    %if (abs(H2) > tol)
    %    error('Far field bc not satisfied')
    %end
    
    v=F1;
    
    figure('position', [0,0,800,800]); 
    plot(eta,v(1,:),'k-','LineWidth',2); hold on; 
    plot(eta,v(2,:),'r-','LineWidth',2); 
    set(gca,'Fontsize',20)
    l1=legend('$v_0(\eta)$','$v_{0\eta}(\eta)$');
    set(l1, 'Interpreter','LaTex','Fontsize',30);
    ylabel('Vel. in the temp. adj. region $v_0$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
    xlim([a,b])
    grid on
    hold off;
    toc
    
end