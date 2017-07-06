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
% [eta, v] = shooting_gotler(@gotler,0.0060,1e-6,1,7,[0 0],'ff',1,0.1);
%
% i.e. solve bvp in gotler on [1,7] with bcs v(1)=0 and v(7)=0, 
% tolerance 1e-6, wavenumbers beta=1 and khat=0.1




%                         Output eigenfuntion                          %



function [eta, v] = shooting_gotler(gotler,deltaeta,tol,a,b,bcs,...
    init,beta,khat) 

    % Parameters and base flow should really be put into funtion 

    gamma=1.4; Pr=1; C=0.509; D=1; etab=1; kappa=0.1;
    
    % Solve for the base flow 
    
    [~,baseT,baseTdash,baseU,intbaseT]= baseflow(C,Pr,D,etab,deltaeta,...
                                        a,b);

    tic; % Begin time
    
    % If my number of arguements is 8 then initial guesses gave been 
    % specified if not take these to be -1 and 1.
    
    if nargin == 10
        shoot1 = init(1); shoot2 = init(2);
    else
        shoot1 = -10; shoot2 = 10;
    end
    
    % Sets up boundary condition vectors, with the first entries being
    % the know dirichlet conditions and the second the two shoots
    
    a1 = [bcs(1) shoot1];
    a2 = [bcs(1) shoot2]; 
    
    % Now iterate solution outwards using Rk method 
    
    [~, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,baseU,kappa,...
        beta,khat,intbaseT); 
    [eta, F2] = RK(a,b,deltaeta,a2,gotler,baseT,baseTdash,baseU,kappa,...
        beta,khat,intbaseT);  
    
    F1 = F1(1,end) - bcs(2);
    F2 = F2(1,end) - bcs(2);
    r = 1;

    % Identify if as root is possible by checking for sign change
    
    % Check 
    % F1*F2
    
    if (F1*F2 > 0) 
        error('The root of F function does not exist')
    end
    
    % Set one shoot for iteration 
    
    F3 = F1;
    
    % Iteration to home in on axis crossing 
    
    while (abs(F3) > tol) 
        
        % Check
        % F3
        
        % Bring one shoot in half the distance between the teo
        
        shoot3 = (shoot1 + shoot2)/2;
        
        % Renforce conditions and rerun RK solver on loop adjusting 
        % to compensate for average overshooting root
        
        a3 = [bcs(1) shoot3];                      
        
        [eta, F3] = RK(a,b,deltaeta,a3,gotler,baseT,baseTdash,baseU,...
            kappa,beta,khat,intbaseT);
        
        % Check
        % F3(r,end);
        
        v = F3; F3 = F3(r,end) - bcs(2); 
        if (F1*F3 < 0)
            shoot2 = shoot3; F2 = F3;            
        elseif (F1*F2 < 0)
            shoot1 = shoot3; F1 = F3;
        else
            error('Something has gone horribly wrong, probs NANS');           
        end
        
    end           
    
    % Plotting of solutions 
    
    % Check 
    % size(eta);
    
    figure('position', [0,0,800,800]); 
    plot(eta,v(1,:),'k-','LineWidth',2); hold on; 
    plot(eta,v(2,:),'r-','LineWidth',2); 
    set(gca,'Fontsize',20)
    l1=legend('$v_0(\eta)$','$v_{0\eta}(\eta)$');
    set(l1, 'Interpreter','LaTex','Fontsize',30);
    ylabel('Vel. in the temp. adj. region $v_0$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
    xlim([1,7])
    grid on
    hold off;
    toc
    
end