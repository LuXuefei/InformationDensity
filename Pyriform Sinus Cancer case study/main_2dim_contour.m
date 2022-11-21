% * Codes for Information Density Estimation using Monte Carlo method
%   The Pyriform Sinus Cancer Case Study of
%  'Information Density in Decision Analysis', 
%   by Gordon B Hazen, Emanuele Borgonovo and Xuefei Lu, 2022
% 
% * Author: Xuefei Lu, xuefei.lu@skema.edu
% * Date: Dec, 2022
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Felli JC, Hazen GB. A Bayesian approach to sensitivity analysis. Health Economics. 1999 May;8(3):263-8.
% [2] Plante DA, Piccirillo JF, Sofferman RA. Decision analysis of treatment options in pyriform sinus carcinoma. Medical Decision Making. 1987 Jun;7(2):74-83.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;close all; clc

inputnames = {'$d_s$','$p_r$','$p_s$','$f_i$','$m_r$','$m_s$', '$q_r$', '$q_s$', '$t_s$'};

PLpara = [0, 0.01, 0.05, 0.1, 1; %1 ds
    0, 0.055, 0.1, 0.15, 1; %2 pr
    0, 0.24, 0.25, 0.55, 1; %3 ps
    0,0,0,0,0;%4 dummy
    0,0,0,0,0;%5 dummy
    0, 0.13, 0.3, 0.5, 1; %6 ms
    0, 0.75, 0.9, 0.98, 1; %7 qr
    0, 0.5, 0.7, 0.95,1]; %8 qs
% Gamma parameters
% multiplier, alpha, beta, power
Gampara = [0.878, 6.392, 1, 0.216;%4 fi
    4.645, 2.041, 1, 0.962];%9 ts
% Uniform parameters
% min, max
Unipara =[0.01, 0.03]; %5 mr

%% base case
ds = 0.05;
pr = 0.1;
ps = 0.25;
fi = 1.3;
mr = 0.01;%%
ms = 0.3;
qr = 0.9;%%
qs = 0.7;
ts = 5;
xbase = [ds, pr, ps, fi, mr, ms, qr, qs, ts];
% 
% [y]=PSC_87(xbase)

%% Find global optimal a^*
Nc = 100000;
X = InputGenerator_Cancer(Nc,PLpara,Gampara,Unipara);
[ya]=PSC_87(X); % use xb for bese senario
% global optimal a^*:
EYa = mean(ya,1); %
[MEYa,astar] = max(EYa);


%% Figure 10 - contour infomation density + optimal area
gidx = [2,7]; % pr qr
xrange = [0,1];
yrange = [0,1];

ngrid = 50;
xgrid = linspace(xrange(1),xrange(2),ngrid);
ygrid = linspace(yrange(1),yrange(2),ngrid);
plotmean = xbase(gidx);
%%
tic
[infogain,astarstar,EXO,EYO ] = EVPI2ddouleloopgrid(xgrid,ygrid,gidx,astar,Nc,PLpara,Gampara,Unipara);
toc 
% ngrid = 50 takes about 36sec
%% joint pdf
Mpdftrue = nan(size(EXO));

for i = 1: ngrid % column
    for j = 1: ngrid  
        x2 = EXO(i,j); x3 = EYO(i,j);%
[Xpdf1] = PSC_87_density(gidx(1),x2,PLpara,Gampara,Unipara);
[Xpdf2] = PSC_87_density(gidx(2),x3,PLpara,Gampara,Unipara);
Mpdftrue(i,j) = Xpdf1*Xpdf2;
    end
end

%%
figure
CO = get(gca,'ColorOrder');
plotdata = infogain.*Mpdftrue;
% contour pdf
c = pink;%hot, bone, pink, white
c = flipud(c);
colormap(c);
contourf(EXO,EYO,plotdata,15)
colorbar
grid on

hold on

[row, col] = find(ismember(plotdata, max(plotdata(:))));
xmax = xgrid(col);
ymax = ygrid(row);
xlabel(inputnames{gidx(1)},'FontSize',20,Interpreter="latex")
ylabel(inputnames{gidx(2)},'FontSize',20,Interpreter="latex")

scatter(xbase(gidx(1)),xbase(gidx(2)),150,'blue','filled','pentagram')
text(xbase(gidx(1))-0.1,xbase(gidx(2))-0.05,['base value'],'HorizontalAlignment','left','Interpreter','latex',FontSize=15,Color='b')
% add optimal
plotdata = astarstar;
contour(EXO,EYO,plotdata,3,'LineWidth',3,'LineColor',CO(2,:))

textposi = [0.7,0.6; 0.3,0.4; 0.05,0.96];%full
text(textposi(1,1),textposi(1,2), ['radiation'],'Interpreter','latex',FontSize=15,Color=CO(2,:))
text(textposi(2,1),textposi(2,2), ['surgery'],'Interpreter','latex',FontSize=15,Color=CO(2,:))
text(textposi(3,1),textposi(3,2), ['surgery then radiation'],'Interpreter','latex',FontSize=15,Color=CO(2,:))

title('$i(p_r,q_r)$',Interpreter='latex',FontSize=15)
