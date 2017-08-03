%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               baseflow                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solves boundary value problem for the base flow using bvp4c. Outputs
% nodes, velocity and temperture fields.

function [eta,baseT,baseTdash,baseU,baseUdash] = ...
          baseflow(C,Pr,D,deltaeta,a,b)

    % matching number
    etab=1;      
    % system of differential equations      
    dydx = @(eta,z)[  z(2)  ;  z(3)  ;  ( -eta*z(3) + ...
        (((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2)) ...
        *z(3) )/(((1+C)*sqrt(z(4)))/(z(4)+C))  ; ...
        z(5)  ;  ( -eta*z(5) + ...
        (Pr^-1)*(((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2))*z(5))...
        /((Pr^-1)*((1+C)*sqrt(z(4)))/(z(4)+C))  ];
    % boundary conditions
    BC = @(za,zb)[  za(1) - D/(etab^(3/Pr))  ;  zb(2)  ; ...
        za(2) + (3/Pr)*D/etab^((3/Pr)-1)  ; ...
        za(4)- (9*(1+C)^2)/(Pr^2*etab^4)  ;  zb(4)-1  ];
    % initial guess   
    zint = @(x)[  0  ;  1  ;  0  ;  1  ;  0 ];
    % initialise solution
    solint=bvpinit(a:deltaeta*10:b,zint);    
    % solve with bvp4c
    S=bvp4c(dydx,BC,solint);
    % Extract solutions for baseflow  
    eta=S.x; baseT=S.y(4,:); baseTdash=S.y(5,:); 
    baseU=S.y(2,:); baseUdash=S.y(3,:);
    % Interpolate base flow for right grid size
    baseT = interp1(eta,baseT,a:deltaeta/5:b,'spline');
    baseTdash = interp1(eta,baseTdash,a:deltaeta/5:b,'spline');
    baseU = interp1(eta,baseU,a:deltaeta/5:b,'spline');
    baseUdash = interp1(eta,baseUdash,a:deltaeta/5:b,'spline');
    eta=a:deltaeta/5:b;

    % plot solution (if running anyother program % out)
    % figure('position', [0,0,800,800]); 
    % plot(eta,baseT,'LineWidth',2); 
    % set(gca,'Fontsize',20)
    % ylabel('Temp. in adj. region, $T_1$','Interpreter',...
    % 'LaTex','Fontsize',40)
    % xlabel('Wall layer variable, $\eta$','Interpreter',...
    % 'LaTex','Fontsize',40)
    % xlim([a,b])
    % grid on 
    % figure('position', [0,0,800,800]); 
    % plot(eta,baseU,'LineWidth',2); 
    % set(gca,'Fontsize',20)
    % ylabel('Vel. in adj. region, $U_1$','Interpreter',...
    % 'LaTex','Fontsize',40)
    % xlabel('Wall layer variable, $\zeta$','Interpreter',...
    % 'LaTex','Fontsize',40)
    % xlim([a,b])
    % grid on
    
end
