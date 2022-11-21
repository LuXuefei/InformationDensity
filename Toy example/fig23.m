% * Codes for Information Density Estimation using Monte Carlo method
%   The Toy Example of
%  'Information Density in Decision Analysis', 
%   by Gordon B Hazen, Emanuele Borgonovo and Xuefei Lu, 2022
% 
% * Author: Xuefei Lu, xuefei.lu@skema.edu
% * Date: Dec, 2022
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clearvars;
%%
pgrid = linspace(0,1,91);
v1 = 10;
v0 = 6;
pcrit = 0.6;
% utility given alterative a
Ua = @(p) max(v1.*p,v0);
% utility improvement
UI = @(p) max(v1.*p,v0) - v1.*p;

% P distribution - beta(a,b)
a=4.8
b=1.2
phat = a/(a+b) % mean 
varp = a*b/((a+b)^2*(a+b+1)) % 0.0229
fpplot = betapdf(pgrid,a,b);

% since v1.*phat > 6, a^* = a1
EUIdp = @(p) (max(v1.*p,v0) - v1.*p).*betapdf(p,a,b)
Infoden = @(p) (max(v1.*p,v0) - v1.*p)
EUIdpplot = EUIdp(pgrid);
Infoplot = Infoden(pgrid);

vert = [pgrid',EUIdpplot'];
f = 1:length(pgrid);

EUI = integral(EUIdp,0,1)
%% Fig 2
figure
CO = get(gca,'ColorOrder');
clf
hold on;
patch('Faces',f,'Vertices',vert,'FaceColor',CO(2,:),'FaceAlpha',.5);
text(0.15,0.75, ['Area = $EUI_p = ', num2str(EUI,2),'$'],'Interpreter','latex',...
    FontSize=20)
ylim([-0.1,2.5])

xline(pcrit,'LineStyle','--',"LineWidth",1)
xline(phat,'LineStyle','-',"LineWidth",2,'Color','b')
[M,I] = max(vert,[],1);
ax = [phat, vert(I(2),1)];
ay = [0.9,0.9]*2.5
set(gcf,'Units','normalized')
set(gca,'Units','normalized')
axx = axis;
ap = get(gca,'Position')
xo = ax;
yo = ay;
xp = (xo-axx(1))/(axx(2)-axx(1))*ap(3)+ap(1);
yp = (yo-axx(3))/(axx(4)-axx(3))*ap(4)+ap(2);
ah=annotation('arrow',xp,yp,'Color','b','LineWidth',2);

xo = [0.25, 0.4];
yo = [0.70, 0.3];
xp = (xo-axx(1))/(axx(2)-axx(1))*ap(3)+ap(1);
yp = (yo-axx(3))/(axx(4)-axx(3))*ap(4)+ap(2);
ae=annotation('arrow',xp,yp);
xlabel('$p$','Interpreter','latex','FontSize',20)
ylabel('$i_{p}(p)$','Interpreter','latex','FontSize',20)

%% Fig 3
% Varies beta parameters while keeping same variance
% https://www.johndcook.com/blog/2021/04/07/beta-given-mean-variance/

pgrid = linspace(0,1,101);
v1 = 10;
v0 = 6;
pcrit = 0.6;

%% Fig 3 - (a)

%a=4.8, b=1.2;
var0 = 0.0229;
mu0 = 0.8;


% figure 4a
mu1 = 0.65;
k = (1-mu1)/mu1;
a = mu1 *(mu1*(1-mu1)/var0 - 1);
b = k*a;

phat = a/(a+b); % mean 
varp = a*b/((a+b)^2*(a+b+1));

% since v1.*phat > 6, a^* = a1
EUIdp = @(p) (max(v1.*p,v0) - v1.*p).*betapdf(p,a,b);
fpplot = betapdf(pgrid,a,b);
EUI = integral(EUIdp,0,1);

%
vert = [pgrid',EUIdp(pgrid)'];
f = 1:length(pgrid);
figure
CO = get(gca,'ColorOrder');
hold on;
patch('Faces',f,'Vertices',vert,'FaceColor',CO(2,:),'FaceAlpha',.5);
title(['$\hat{p} = ',num2str(mu1),', EUI_p =', num2str(EUI,2),'$'],'Interpreter','latex','FontSize',18)
ylim([-0.1,2.5])

xline(pcrit,'LineStyle','--',"LineWidth",1)
xline(phat,'LineStyle','-',"LineWidth",2,'Color','b')
[M,I] = max(vert,[],1);
ax = [phat, vert(I(2),1)];
ay = [0.9,0.9]*2.5
set(gcf,'Units','normalized')
set(gca,'Units','normalized')
axx = axis;
ap = get(gca,'Position')
xo = ax;
yo = ay;
xp = (xo-axx(1))/(axx(2)-axx(1))*ap(3)+ap(1);
yp = (yo-axx(3))/(axx(4)-axx(3))*ap(4)+ap(2);
ah=annotation('arrow',xp,yp,'Color','b','LineWidth',2);
xlabel('$p$','Interpreter','latex','FontSize',20)
ylabel('$i_{p}(p)$','Interpreter','latex','FontSize',20)
%% Fig 3 - (b)
mu1 = 0.55;
k = (1-mu1)/mu1;
a = mu1 *(mu1*(1-mu1)/var0 - 1);
b = k*a;

phat = a/(a+b); % mean 
varp = a*b/((a+b)^2*(a+b+1));

% since v1.*phat < v0, a^* = a0
EUIdp = @(p) (max(v1.*p,v0) - v0).*betapdf(p,a,b);
fpplot = betapdf(pgrid,a,b);
EUI = integral(EUIdp,0,1)

%
vert = [pgrid',EUIdp(pgrid)'];
f = 1:length(pgrid);
figure
CO = get(gca,'ColorOrder');
hold on;
patch('Faces',f,'Vertices',vert,'FaceColor',CO(2,:),'FaceAlpha',.5);
title(['$\hat{p} = ',num2str(mu1),', EUI_p =', num2str(EUI,2),'$'],'Interpreter','latex','FontSize',18)
ylim([-0.1,2.5])


xline(pcrit,'LineStyle','--',"LineWidth",1)
xline(phat,'LineStyle','-',"LineWidth",2,'Color','b')
[M,I] = max(vert,[],1);
ax = [phat+0.01, vert(I(2),1)];
ay = [0.9,0.9]*2.5
set(gcf,'Units','normalized')
set(gca,'Units','normalized')
axx = axis;
ap = get(gca,'Position')
xo = ax;
yo = ay;
xp = (xo-axx(1))/(axx(2)-axx(1))*ap(3)+ap(1);
yp = (yo-axx(3))/(axx(4)-axx(3))*ap(4)+ap(2);
ah=annotation('arrow',xp,yp,'Color','b','LineWidth',2);
xlabel('$p$','Interpreter','latex','FontSize',20)
ylabel('$i_{p}(p)$','Interpreter','latex','FontSize',20)

%% Fig 3 - (c)
mu1 = 0.40;
k = (1-mu1)/mu1;
a = mu1 *(mu1*(1-mu1)/var0 - 1);
b = k*a;

phat = a/(a+b); % mean 
varp = a*b/((a+b)^2*(a+b+1));

% since v1.*phat < v0, a^* = a0
EUIdp = @(p) (max(v1.*p,v0) - v0).*betapdf(p,a,b);
fpplot = betapdf(pgrid,a,b);
EUI = integral(EUIdp,0,1)

%
vert = [pgrid',EUIdp(pgrid)'];
f = 1:length(pgrid);
figure
CO = get(gca,'ColorOrder');
hold on;
patch('Faces',f,'Vertices',vert,'FaceColor',CO(2,:),'FaceAlpha',.5);

title(['$\hat{p} = ',num2str(mu1),', EUI_p =', num2str(EUI,2),'$'],'Interpreter','latex','FontSize',18)
ylim([-0.1,2.5])


xline(pcrit,'LineStyle','--',"LineWidth",1)
xline(phat,'LineStyle','-',"LineWidth",2,'Color','b')
[M,I] = max(vert,[],1);
ax = [phat+0.02, vert(I(2),1)];
ay = [0.9,0.9]*2.5
set(gcf,'Units','normalized')
set(gca,'Units','normalized')
axx = axis;
ap = get(gca,'Position')
xo = ax;
yo = ay;
xp = (xo-axx(1))/(axx(2)-axx(1))*ap(3)+ap(1);
yp = (yo-axx(3))/(axx(4)-axx(3))*ap(4)+ap(2);
ah=annotation('arrow',xp,yp,'Color','b','LineWidth',2);
xlabel('$p$','Interpreter','latex','FontSize',20)
ylabel('$i_{p}(p)$','Interpreter','latex','FontSize',20)

