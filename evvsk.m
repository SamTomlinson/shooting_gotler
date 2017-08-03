%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 evvsk                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runs shooting_gotler3 for a variety of different wavenumbers then
% plots the dependance of the ev on k and some of the different 
% eigenmodes

%% calculate wavenumber dependancy 

% flow parameters 
a=0.5; b=20; deltaeta=0.01;
Pr=1; C=0.509; D=1; A=3*(1+C)/Pr;
% base flow
[~,baseT,baseTdash,baseU,baseUdash]=baseflow(C,Pr,D,deltaeta,a,b);
%initialise 
ev=[];
% loop through different wavenumbers and calculate ev
for k=0.1:0.1:8
    k
    % calculate ev
    [eta, v,eigval] = shooting_gotler3(@gotler,deltaeta,a,b,k);
    % fill vector
    ev=[eigval, ev];
    % extract ones to plot
    if k==0.5
        v1=v(1,:);
        dv1=v(2,:);
    end
    if k==1
        v2=v(1,:);
        dv2=v(2,:);
    end
    if k==2
        v3=v(1,:);
        dv3=v(2,:);
    end     
end
% reverse ev order as come out the wrong way round
ev=ev(end:-1:1);

%% normalisation of other flow variables

% recalculate base eta and base flow variables with correct size
baseeta=a:deltaeta/5:b;
baseTm = interp1(baseeta,baseT,a:deltaeta:b,'spline');
baseTdashm = interp1(baseeta,baseTdash,a:deltaeta:b,'spline');
baseUm = interp1(baseeta,baseU,a:deltaeta:b,'spline');
baseUdashm = interp1(baseeta,baseUdash,a:deltaeta:b,'spline');
T1=normalise(-baseTdashm.*v1./baseTm);
T2=normalise(-baseTdashm.*v2./baseTm); 
T3=normalise(-baseTdashm.*v3./baseTm);
u1=-baseUdashm.*v1./baseTm+1; 
u2=-baseUdashm.*v2./baseTm+1; 
u3=-baseUdashm.*v3./baseTm+1;
w1=normalise(-dv1./(0.5.*baseTm)); 
w2=normalise(-dv2./(1.*baseTm)); 
w3=normalise(-dv3./(2.*baseTm));
p1=normalise(-dv1./(0.5.^2.*baseTm.^2));
p2=normalise(-dv2./(1.^2.*baseTm.^2)); 
p3=normalise(-dv3./(2.^2.*baseTm.^2));

%% plotting

figure('position', [0,0,800,800]); 
plot(baseeta,baseT,'LineWidth',2); 
set(gca,'Fontsize',20)
ylabel('Temp. in adj. region, $T_1$','Interpreter', 'LaTex','Fontsize',40)
xlabel('Wall layer variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
xlim([0,b])
grid on
figure('position', [0,0,800,800]); 
plot(baseeta,baseU,'LineWidth',2); 
set(gca,'Fontsize',20)
ylabel('Vel. in adj. region, $U_1$','Interpreter', 'LaTex','Fontsize',40)
xlabel('Wall layer variable, $\zeta$','Interpreter', 'LaTex','Fontsize',40)
xlim([0,b])
grid on
figure('position', [0,0,800,800]); 
plot(0.1:0.5:8,ev); hold on;
set(gca,'Fontsize',20)
ylabel('Eigenvalue, $\beta^2(G^*-Q)$','Interpreter',...
        'LaTex','Fontsize',40)
xlabel('Wavenumber, $\hat{k}$','Interpreter', 'LaTex','Fontsize',40)
xlim([0.5,8])
grid on
hold off;
figure('position', [0,0,800,800]); 
plot(eta,v1,'LineWidth',2); hold on; 
plot(eta,v2,'LineWidth',2);  
plot(eta,v3,'LineWidth',2); 
set(gca,'Fontsize',20)
l1=legend('$k=0.5$','$k=1$','$k=2$');
set(l1, 'Interpreter','LaTex','Fontsize',30);
ylabel('Vel. in the temp. adj. region $v_0$','Interpreter',...
        'LaTex','Fontsize',40)
xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
xlim([0,b])
grid on
hold off;







