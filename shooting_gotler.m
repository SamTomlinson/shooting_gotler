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




%                               Example                                %
%
% [eta, v] = shooting_gotler(@gotler,deltaeta,tol,a,b,bcs,init,beta,khat)
% [eta, v] = shooting_gotler(@gotler,0.0060,1e-6,1,7,[0 0],1,0.1);
%
% i.e. solve bvp in gotler on [1,7] with bcs v(1)=0 and v(7)=0, 
% tolerance 1e-6, wavenumbers beta=1 and khat=0.1




%                         Output eigenfuntion                          %



function [eta, v] = shooting_gotler(gotler,deltaeta,tol,a,b,eigval,init) 

    % Parameters and base flow should really be put into funtion 

    gamma=1.4; Pr=1; C=0.509; D=1; etab=1;
    eta=a:deltaeta:b; A=3*(1+C)/Pr;
    
    
    % Solve for the base flow 
    
    [~,baseT,baseTdash]=baseflow(C,Pr,D,1,deltaeta,a,b);

    tic; % Begin time
    
    % Initial shoots
    
    shoot1 = 0.1; shoot2 = 8;
    
    % Sets up far field conditions
    
    a1 = [exp(-shoot1*b), -shoot1*exp(shoot1*b)];
    a2 = [exp(-shoot2*b), -shoot2*exp(shoot2*b)];
    
    
    % Now iterate solution backwards using Rk method 
    
    [~, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,shoot1,eigval); 
    [eta, F2] = RK(a,b,deltaeta,a2,gotler,baseT,baseTdash,shoot2,eigval);
    
    F1;
    
    % Calculate boundary decay 
    
    H1=F1(2,1)-(shoot1*A^2/(a^4))*F1(1,1);
    H2=F2(2,1)-(shoot2*A^2/(a^4))*F2(1,1);
    
    % Notify of errors
    
    if (H1*H2 > 0) 
        error('The root does not exist')
    end
    
    % Set one shoot for iteration 
    
    H3 = H1;
    
    % Iteration to home in on axis crossing
    
%     while (abs(H3) > tol) 
%         
%         H1;
%         H2;
%        
         shoot3=(shoot1+shoot2)/2;
% 
%         Renforce conditions and rerun RK solver on loop adjusting 
%         to compensate for average overshooting root
%         
         a3 = [exp(-shoot3*b), -shoot3*exp(shoot3*b)];                     
%         
         [eta, F3] = RK(a,b,deltaeta,a3,gotler,baseT,baseTdash,shoot3,eigval);
%         
%         Check
%         F3(r,end);
%         
         v = F3; 
%         
%         H3=F3(2,1)-(shoot3*A^2/(a^4))*F3(1,1)
%     
%         if (H3 < 0)
%             shoot1 = shoot3;
%             F1 = F3;  
%             H1 = H3;
%         elseif (H3 > 0)
%             shoot2 = shoot3;
%             F2 = F3;
%             H2 = H3;
%         else
%             error('Something has gone horribly wrong, probs NANS');           
%         end
%         
%         size(v);
%         size(eta);
%         
%         khat=shoot3;
%         
%     end 
    
    % Different sort of attempt
    
    vec=[];
    kvec=[];
    for shoot1=-1:0.001:1
    a1 = [exp(-shoot1*b), -shoot1*exp(-shoot1*b)];
    [~, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,shoot1,eigval);     
    H1=F1(2,1)-((shoot1*A^2)/(a^4))*F1(1,1);
    H2=F1(2,length(F1(2,:))) + shoot1*F1(1,length(F1(2,:)));
    vec=[H1,vec];
    kvec=[shoot1,kvec];
    if shoot1==-0.5
    display('25percentage completed')
    elseif shoot1==0
    display('50percentage completed')
    elseif shoot1==0.5
    display('75percentage completed')
    elseif shoot1==1
    display('100percentage completed')
    end 
    end
    plot(kvec,vec);
    tmp = abs(vec-0)
    [idx, idx] = min(tmp) %index of closest value
    closest = vec(idx)
    kval= kvec(idx)
    
    % Plotting of solutions 
    
    % Check 
    % size(eta);
    
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
    toc
    
end