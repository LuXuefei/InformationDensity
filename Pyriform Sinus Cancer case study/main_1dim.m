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

clearvars;close all
clc

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

% For a more precised estimation, use Nc = 100000.
Nc = 10000;
X = InputGenerator_Cancer(Nc,PLpara,Gampara,Unipara);

%% Base case simulation
ds = 0.05;
pr = 0.1;
ps = 0.25;
fi = 1.3;
mr = 0.01;%%
ms = 0.3;
qr = 0.9;%%
qs = 0.7;
ts = 5;
xbasevalues = [ds, pr, ps, fi, mr, ms, qr, qs, ts];

[y]=PSC_87(xbasevalues) % 72.4261   97.1935  103.1286  108.1899

%% Generate dataset [x,ya]
[ya]=PSC_87(X); % use xb for bese senario
% global optimal a^*:
EYa = mean(ya,1) % use xb for bese senario, EYa = Table5.
[MEYa,astar] = max(EYa)

% total VOI
ET = mean(max(ya,[],2)) - mean(ya(:,astar)) % EVPI/ET --> for display

%% Double-loop Monte-Carlo estimation EVPI - Table 1. EUI
% When use Nc=10000, it takes about 206 sec
% When use Nc = 100000, it takes about 6220sec
tic
EVPId = EVPIdouleloop(X)
toc

figure
bar(1:9,EVPId);
for i1=1:numel(EVPId)
    text(i1,EVPId(i1),num2str(EVPId(i1),'%0.4f'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',14)
end
title('EUI')
set(gca,'XTickLabel',inputnames,'FontSize',13,'TickLabelInterpreter','latex');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Double loop: 1-dim alternative/infomation gain/info density
% Figures 7 - 9
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the input index
i = 4 % From 1 to 9

ngrid = 100;

if   ismember(i,[1:3,6:8])
    xgrid = linspace(PLpara(i,1),PLpara(i,5),ngrid);
    plotmean = PLpara(i,3); % basevalue
elseif i==4
    xgrid = linspace(0,1.8,ngrid);
    plotmean = 1.3; % basevalue
elseif i == 9
    xgrid = linspace(0,50,ngrid);
    plotmean = 5; % basevalue
elseif i== 5
    xgrid = linspace(Unipara(1),Unipara(2),ngrid);
    plotmean = 0.01; % basevalue
end

astarstar = nan(size(xgrid)); % optimal alternative
EYAstarstarXu = nan(size(xgrid));
EYAstarXu = nan(size(xgrid));


for k=1:length(xgrid)
    xu = X;
    xu(:,i) = xgrid(k);
    [yaxu]=PSC_87(xu);
    EYaXu = mean(yaxu);
    [EYAstarstarXu(k),astarstar(k)] = max(EYaXu);
    EYAstarXu(k) = EYaXu(astar);
end
infogain = EYAstarstarXu - EYAstarXu;

figure
CO = get(gca,'ColorOrder');
subplot(2,1,1)
yyaxis left
plot(xgrid,astarstar,'Color',CO(5,:),'LineWidth',2);
ylim([0.999,4.001])
ylabel('$a^{**}(x)$',Interpreter="latex")

yyaxis right
plot(xgrid,infogain,'Color',CO(2,:),'LineWidth',2); hold on
ylabel(['$\zeta$(x)'],Interpreter="latex")
xline(plotmean,'Color','b','LineWidth',1.5);
if ~isempty(xgrid(diff(astarstar)~=0))
    xline(xgrid(diff(astarstar)~=0),'--')
end
xlabel(inputnames{i},Interpreter="latex")
%legend({'Opt alternative','Info gain'},Location="best")
title('Optimal alternative, Infomation gain')

subplot(2,1,2)
[Xpdf] = PSC_87_density(i,xgrid,PLpara,Gampara,Unipara);
infoden = infogain.*Xpdf;
plot(xgrid,infoden,'Color',CO(1,:),'LineWidth',2)
xlabel(inputnames{i},Interpreter="latex")
title('Infomation density')
xline(plotmean,'Color','b','LineWidth',1.5);
if ~isempty(xgrid(diff(astarstar)~=0))
    xline(xgrid(diff(astarstar)~=0),'--')
end
ylabel('$i_(x)$',Interpreter='latex')

set(gcf, 'PaperPosition', [0 0 7 8]/2); %[left bottom width height]
set(gcf, 'PaperSize', [7 8]/2);
if ~ismember(i,[5,6,9])
    [~,ii] = max(infoden);
    if xgrid(ii)>plotmean
        axend = plotmean + 0.5*(xgrid(ii)-plotmean);
    else
        axend = plotmean -0.5* (plotmean-xgrid(ii));
    end
    ax = [plotmean, axend];
    ay = 0.9*[max(get(gca,'Ylim')),max(get(gca,'Ylim'))];
    set(gcf,'Units','normalized')
    set(gca,'Units','normalized')
    axx = axis;
    ap = get(gca,'Position');
    xo = ax;
    yo = ay;
    xp = (xo-axx(1))/(axx(2)-axx(1))*ap(3)+ap(1);
    yp = (yo-axx(3))/(axx(4)-axx(3))*ap(4)+ap(2);
    ah=annotation('arrow',xp,yp,'Color','b','LineWidth',2);
end

