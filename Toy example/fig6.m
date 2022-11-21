% * Codes for Information Density Estimation using Monte Carlo method
%   The Toy Example of
%  'Information Density in Decision Analysis', 
%   by Gordon B Hazen, Emanuele Borgonovo and Xuefei Lu, 2022
% 
% * Author: Xuefei Lu, xuefei.lu@skema.edu
% * Date: Dec, 2022
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;close all; clc

%%
v0=65;
% P distribution - beta(a,b)
a=4.8;
b=1.2;
phat = a/(a+b); % mean 
varp = a*b/((a+b)^2*(a+b+1));
pgrid = linspace(0,1,91);
fpplot = betapdf(pgrid,a,b);

% V1 distribtuion - normal (80, 10)
mu_v = 80;
sigma_v = 10;
vgrid = linspace(-50,200,101);
fvplot = normpdf(vgrid, mu_v,sigma_v);
%% Fig5 (a) - 2d pdf

figure
[X,Y] = meshgrid(pgrid,vgrid);
Z = betapdf(X,a,b).*normpdf(Y, mu_v,sigma_v);

% contour pdf
c = white;%hot, bone, pink, white
c = flipud(c);
colormap(c);
contourf(X,Y,Z,15)


DecBoundaryV1 = @(p) v0./p;
boundarycolor = 0.2*[1,1,1]; %CO(2,:);
hold on
plot(pgrid,DecBoundaryV1(pgrid),'--',"LineWidth",2,'Color',boundarycolor)
ylim([-40,200])
text(0.03,150,'$p\times v_1 = \$ 65$','Interpreter','latex',FontSize=15,Color=boundarycolor)
text(0.7,190,'$a^*(p,v_1) = a_1$','Interpreter','latex',FontSize=15,Color=CO(2,:))
text(0.2,50,'$a^*(p,v_1) = a_0$','Interpreter','latex',FontSize=15,Color=CO(2,:))
xlabel('$p$','Interpreter','latex','FontSize',20)
ylabel('$v_1$','Interpreter','latex','FontSize',20)
title('Joint pdf of $p$ and $v_1$','Interpreter','latex','FontSize',20)
grid on

%% Fig5 b- 2d information density

figure
[X,Y] = meshgrid(pgrid,vgrid);
Z = betapdf(X,a,b).*normpdf(Y, mu_v,sigma_v);
% since v1.*phat = 64 < v0 = 65, a^* = a0
EUIdpv1plot = (max(X.*Y,v0) - v0).*Z

EUIdpv1 = @(p,v) (max(p.*v,v0) - v0).* betapdf(p,a,b).*normpdf(v, mu_v,sigma_v)
EUI = integral2(EUIdpv1,0,1,-40,200);

% contour pdf
c = pink;%hot, bone, pink, white
c = flipud(c);
colormap(c);
contourf(X,Y,EUIdpv1plot,15)
[row, col] = find(ismember(EUIdpv1plot, max(EUIdpv1plot(:))))
vmax = vgrid(row);
pmax = pgrid(col);
%
DecBoundaryV1 = @(p) v0./p;
boundarycolor = 0.2*[1,1,1]; %CO(2,:);
hold on
plot(pgrid,DecBoundaryV1(pgrid),'--',"LineWidth",2,'Color',boundarycolor)
ylim([-40,200])

xlabel('$p$','Interpreter','latex','FontSize',20)
ylabel('$v_1$','Interpreter','latex','FontSize',20)
title('$i_{p,v_1}$','Interpreter','latex','FontSize',20)
grid on
ax = [0.80, pmax];
ay = [80,vmax];
set(gcf,'Units','normalized')
set(gca,'Units','normalized')
axx = axis;
ap = get(gca,'Position')
xo = ax;
yo = ay;
xp = (xo-axx(1))/(axx(2)-axx(1))*ap(3)+ap(1);
yp = (yo-axx(3))/(axx(4)-axx(3))*ap(4)+ap(2);
ah=annotation('arrow',xp,yp,'Color','b','LineWidth',2);

text(0.8,120,['(',num2str(pmax,2),',',num2str(vmax,3),')'],'Interpreter','latex',FontSize=15,Color='b')
text(0.5,80,['(',num2str(0.8),',',num2str(80),')'],'Interpreter','latex',FontSize=15,Color='b')


