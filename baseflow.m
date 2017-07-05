function [x,baseT,baseTdash,baseU,intbaseT]= baseflow(C,Pr,D,eta)
    
dydx=@(x,z)[z(2);z(3);(-x*z(3)+ (((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2))*z(3))/(((1+C)*sqrt(z(4)))/(z(4)+C)); ...
    z(5);(-x*z(5) + (Pr^-1)*(((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2))*z(5))/((Pr^-1)*((1+C)*sqrt(z(4)))/(z(4)+C))];

BC=@(za,zb)[za(1) - D/(eta^(3/Pr)) ; zb(2) ; za(2) + (3/Pr)*D/eta^((3/Pr)-1); za(4)- (9*(1+C)^2)/(Pr^2*eta^4); zb(4)-1];
    
zint=@(x)[0 ; 1; 0 ; 1 ; 0];
    
solint=bvpinit(linspace(1,3,1001),zint);
    
S=bvp4c(dydx,BC,solint);
    
x=S.x; baseT=S.y(4,:); baseTdash=S.y(5,:); intbaseT=trapz(baseT);
 
baseU=S.y(2,:);

figure('position', [0,0,800,800]); 
plot(x,baseT,'LineWidth',2); 
set(gca,'Fontsize',20)
%hold on; plot(a2,b2); plot(a3,b3); plot(a4,b4); plot(a5,b5); 
%legend('n=0.2 Num','n=0.4  Num','n=0.6  Num','n=0.8  Num','n=1.0  Num')
ylabel('Temp. in adj. region, $T_1$','Interpreter', 'LaTex','Fontsize',40)
xlabel('Wall layer variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
xlim([1,3])
grid on

figure('position', [0,0,800,800]); 
plot(x,baseU,'LineWidth',2); 
set(gca,'Fontsize',20)
%hold on; plot(a2,c2); plot(a3,c3); plot(a4,c4); plot(a5,c5)
%legend('n=0.2 Num','n=0.4  Num','n=0.6  Num','n=0.8  Num','n=1.0  Num')
ylabel('Vel. in adj. region, $U_1$','Interpreter', 'LaTex','Fontsize',40)
xlabel('Wall layer variable, $\zeta$','Interpreter', 'LaTex','Fontsize',40)
xlim([1,3])
grid on
    