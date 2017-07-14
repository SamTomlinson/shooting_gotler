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



function [eta, v,eigval] = shooting_gotler3(gotler,deltaeta,a,b,khat) 

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
    
%     figure('position', [0,0,800,800]); 
%     plot(eigvec,vec,'k-','LineWidth',2); 
%     set(gca,'Fontsize',20)
%     ylabel('Near field error, $H(\hat{k})$','Interpreter',...
%         'LaTex','Fontsize',40)
%     xlabel('Eig. val., $\beta^2(G^*-Q)$','Interpreter', 'LaTex','Fontsize',40)
%     xlim([0.01,0.6])
%     %ylim([-1,1])
%     grid on
%     hold off;
   
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
% eigvalright=eigvec(zerIdx(1)-1);
% eigvalleft=eigvec(zerIdx(1)+1);
% vecsright=vec(zerIdx(1)-1);
% vecsleft=vec(zerIdx(1)+1);
eigval=eigs(1);


% eigvalright=0.342; eigvalleft=0.34;  vecsright=-1.426036221091628e+08;
% vecsleft=1.677088912688837e+06;

% Iterate in 

% while (vecsleft>1e-16)
% eignew=(eigvalleft+eigvalright)/2;
% vecnew=(vecsleft+vecsright)/2;
% if vecnew<0
%     eigvalright=eignew;
%     vecsright=vecnew;
% end
% if vecnew>0
%     eigvalleft=eignew;
%     vecsleft=vecnew;
% end
% end
% 
% eigval=eigvalleft;

% Improve accuracy 
diff=1;

tol=0.1;
while abs(diff>1e-16)
    eigvalold=eigval;
    [eigval,H1]=loop(eigval,khat,a,b,A,deltaeta,a1,gotler,baseT,...
    baseTdash,shoot1,tol);
    diff=abs(eigvalold-eigval);
    tol=tol/10;
end


       
% Calculation of eigenmodes 
    
a1 = [exp(-khat*b), -khat*exp(-khat*b)];
[eta, F1] = RK(a,b,deltaeta,a1,gotler,baseT,baseTdash,khat,eigval);     
H1=F1(2,1)-((khat*A^2)/(a^4))*F1(1,1);
H2=F1(2,length(F1(2,:))) + khat*F1(1,length(F1(2,:)));
v=F1;

% Set all values below zero to zero (not saure if this is allowed_

v1=v(1,:);
v1(v1<0)=0;
v(1,:)=v1;

% Find local maxima

maxs=[];
for j=2:length(v(1,:))-1
    if (v(1,j-1)<v(1,j) && v(1,j)>v(1,j+1))
        maxs=[j,maxs];
    end
end

% Normalise

v(1,:)=(1/v(1,maxs(1)))*v(1,:);


% Normalisation of eigenmodes



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
%     
%     figure('position', [0,0,800,800]); 
%     plot(eta,-baseTdash.*v(1,:)./baseT,'k-','LineWidth',2);
%     set(gca,'Fontsize',20)
%      ylabel('Temp. in the temp. adj. region $T_0$','Interpreter',...
%         'LaTex','Fontsize',40)
%     xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
%     xlim([0.1,b])
%     grid on
%     hold off;
    toc
    
end