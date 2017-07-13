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
% eigval - found eigenvalue for entered khat
%
% k - spanwise wavenumber
% 
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE



%                         Output eigenfuntion                          %



function [eta, v,eignew1] = shooting_gotler3(gotler,deltaeta,a,b,khat) 

    % Parameters and base flow should really be put into funtion 

    Pr=1; C=0.509; D=1; A=3*(1+C)/Pr; tol=1e-6;
    
    % Solve for the base flow 
    
    [~,baseT,baseTdash,~,~]=baseflow(C,Pr,D,deltaeta,a,b);

    % Begin time
    
    tic;
    
    % Initialise vectors
    
    vec=[]; eigvec=[];
    
    % Loop through different khat values 
    
    for shoot1=0.01:0.001:0.6
        
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
    
    
% Plot H vs eig
    
    figure('position', [0,0,800,800]); 
    plot(eigvec,vec,'k-','LineWidth',2); 
    set(gca,'Fontsize',20)
    ylabel('Near field error, $H(\hat{k})$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('Eig. val., $\beta^2(G^*-Q)$','Interpreter', 'LaTex','Fontsize',40)
    xlim([0.01,0.6])
    %ylim([-1,1])
    grid on
    hold off;
   
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
grad=(vecsleft-vecsright)/(eigvalleft-eigvalright);
angle=pi/2+atan(grad);
radtodeg(angle);
deltaev=vecsleft*tan(angle);
delta2ev=vecsright*tan(angle);
eignew1=eigvalleft+deltaev;
eignew2=eigvalright+delta2ev;
   
% Calculation of eigenmodes 
    
a1 = [exp(-khat*b), -khat*exp(-khat*b)];
[eta, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,khat,eignew1);     
H1=F1(2,1)-((khat*A^2)/(a^4))*F1(1,1)
H2=F1(2,length(F1(2,:))) + khat*F1(1,length(F1(2,:)))
v=F1;

% Plotting of eigenomdes (if running evvsk % out)

    figure('position', [0,0,800,800]); 
    plot(eta,v(1,:),'LineWidth',2); 
    set(gca,'Fontsize',20)
    ylabel('Vel. in the temp. adj. region $v_0$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
    xlim([0.1,b])
    grid on
    hold off;
    
    figure('position', [0,0,800,800]); 
    plot(eta,-baseTdash.*v(1,:)./baseT,'k-','LineWidth',2);
    set(gca,'Fontsize',20)
     ylabel('Temp. in the temp. adj. region $T_0$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
    xlim([0.1,b])
    grid on
    hold off;
    toc
    
end