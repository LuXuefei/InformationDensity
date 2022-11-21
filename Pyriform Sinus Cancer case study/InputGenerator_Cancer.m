function X = InputGenerator_Cancer(Nc,PLpara,Gampara,Unipara)
% Table 1 Felli JC, Hazen GB. A Bayesian approach to sensitivity analysis. Health Economics. 1999 May;8(3):263-8.

plotfigure = 0; % or 0
inputnames = {'$d_s$','$p_r$','$p_s$','$f_i$','$m_r$','$m_s$', '$q_r$', '$q_s$', '$t_s$'};
% %Piecewise linear parameters
% PLpara = [0, 0.01, 0.05, 0.1, 1; %1 ds
%     0, 0.055, 0.1, 0.15, 1; %2 pr
%     0, 0.24, 0.25, 0.55, 1; %3 ps
%     0,0,0,0,0;%4 dummy
%     0,0,0,0,0;%5 dummy
%     0, 0.13, 0.3, 0.5, 1; %6 ms
%     0, 0.75, 0.9, 0.98, 1; %7 qr
%     0, 0.5, 0.7, 0.95,1]; %8 qs
% % Gamma parameters
% % multiplier, alpha, beta, power
% Gampara = [0.878, 6.392, 1, 0.216;%4 fi
%     4.645, 2.041, 1, 0.962];%9 ts
% % Uniform parameters
% % min, max
% Unipara =[0.01, 0.03]; %5 mr


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

X = zeros(Nc,9);


%Calculate piecewise linear pdf parameters
fa1 = @(l,n) 0.025*2/((l-n)*l);
fa4 = @(m,u) -2*0.025/((m-u)^2);
fb4 = @(m,a4) -m*a4;
fa2 = @(b,u,l,m,a1,a4) 1/((u-l)*(b-l)) * (0.95*2 - a4*(u-b)*(u-m) - a1*l*(b-l)) - a1*l/(b-l);
fb2 = @(a1,a2,l) (a1-a2)*l;
fa3 = @(u,b,a4,b4,a2,b2) 1/(u-b)*(a4*u + b4 - b2 - a2*b);
fb3 = @(a2,b,b2,a3) a2*b + b2 - a3*b;

%%
if plotfigure == 1
    figure
    CO = get(gca,'ColorOrder');
end


for idx = [1:3,6:8]
    xmin = PLpara(idx,1);
    xlow = PLpara(idx,2);
    xbase = PLpara(idx,3);
    xupp = PLpara(idx,4);
    xmax = PLpara(idx,5);

    % must be in the following order
    a1 = fa1(xlow,xmin);
    a4 = fa4(xmax,xupp);
    b4 = fb4(xmax,a4);
    % for x8 qr, triangle pdf
    if (a4 == Inf) || (b4 == Inf)
        a4=0;b4=0;
        a1 = a1*2; % cumulative prob [min,low]  = 0.025*2 = 0.05
    end
    a2 = fa2(xbase,xupp,xlow,xmax,a1,a4);
    b2 = fb2(a1,a2,xlow);
    a3 = fa3(xupp,xbase,a4,b4,a2,b2);
    b3 = fb3(a2,xbase,b2,a3);

    fpiecewise = @(x) (a1*x).*(xmin<= x & x<xlow)+ (a2*x+ b2).*(xlow<= x & x<xbase) + ...
        (a3*x+ b3).*(xbase<= x & x<xupp) + (a4*x+ b4).*(xupp<= x & x<xmax);

    X(:,idx) = sampleDist(fpiecewise, 20,Nc,[xmin,xmax]); % func, threshold, n, [xmin,xmax], plot

    if plotfigure == 1
        % test and plot
        subplot(3,4,idx)
        hold on
        %histogram(X(:,idx),Normalization="pdf");
        xline(mean(X(:,idx)),'Color','b','LineWidth',1.5,'LineStyle','--');
        xgrid = linspace(xmin,xmax,100);
        yplot = fpiecewise(xgrid);
        plot(xgrid,yplot, 'LineWidth',2, 'Color',CO(2,:))
        xline(xbasevalues(idx),'Color','b','LineWidth',1.5);
        title(inputnames{idx},Interpreter="latex",FontSize=15)
        ylabel('pdf')
        testpdf=integral(fpiecewise,xmin,xmax);
    end

end
%% Gamma
% X5-fi
idx = 4;
X(:,idx) = Gampara(1,1).*gamrnd(Gampara(1,2),Gampara(1,3),Nc,1).^Gampara(1,4);

a = Gampara(1,1);
alpha = Gampara(1,2);
beta = Gampara(1,3);
b = Gampara(1,4);
fgamma = @(y) gampdf((y./a).^(1/b),alpha,beta).* (1/(a*b)).*(y./a).^(1/b-1);
%testpdf=integral(fgamma,0,1.8);

if plotfigure == 1
    subplot(3,4,idx)
    hold on
    %histogram(X(:,idx),Normalization="pdf");
    xline(mean(X(:,idx)),'Color','b','LineWidth',1.5,'LineStyle','--');
    xgrid = linspace(0,1.8,100);
    yplot = fgamma(xgrid);
    plot(xgrid,yplot, 'LineWidth',2, 'Color',CO(2,:))
    xline(xbasevalues(idx),'Color','b','LineWidth',1.5);
    title(inputnames{idx},Interpreter="latex",FontSize=15)
    ylabel('pdf')
end


% X10, ts
idx = 9;
X(:,idx) = Gampara(2,1).*gamrnd(Gampara(2,2),Gampara(2,3),Nc,1).^Gampara(2,4);

a = Gampara(2,1);
alpha = Gampara(2,2);
beta = Gampara(2,3);
b = Gampara(2,4);

fgamma = @(y) gampdf((y./a).^(1/b),alpha,beta).* (1/(a*b)).*(y./a).^(1/b-1);

if plotfigure == 1
    subplot(3,4,idx)
    hold on
    %histogram(X(:,idx),Normalization="pdf");
    xline(mean(X(:,idx)),'Color','b','LineWidth',1.5,'LineStyle','--');
    xgrid = linspace(0,50,100);
    yplot = fgamma(xgrid);
    plot(xgrid,yplot, 'LineWidth',2, 'Color',CO(2,:))
    xline(xbasevalues(idx),'Color','b','LineWidth',1.5);
    
    title(inputnames{idx},Interpreter="latex",FontSize=15)
    ylabel('pdf')
end

%% Uniform

% X6 - mr
idx = 5;
X(:,idx) = Unipara(1) + (Unipara(2) -Unipara(1) )*rand(Nc,1);

if plotfigure == 1
    subplot(3,4,idx)
    hold on
    %histogram(X(:,idx),Normalization="pdf");
    xline(mean(X(:,idx)),'Color','b','LineWidth',1.5,'LineStyle','--');
    xgrid = linspace(Unipara(1),Unipara(2),100);
    yplot = 1/(Unipara(2) - Unipara(1))*ones(size(xgrid));
    plot(xgrid,yplot, 'LineWidth',2, 'Color',CO(2,:))
    xline(xbasevalues(idx),'Color','b','LineWidth',1.5);
    title(inputnames{idx},Interpreter="latex",FontSize=15)
    ylabel('pdf')

    set(gcf, 'PaperPosition', [-1 -0.5 12 7]); %[left bottom width height]
    set(gcf, 'PaperSize', [10 6.5]);
    saveas(gcf, ['Cancer_inputhistdis_pdf.pdf'])
end



end
