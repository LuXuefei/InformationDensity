function [Xpdf] = PSC_87_density(gidx,xgrid,PLpara,Gampara,Unipara)

inputnames = {'$d_s$','$p_r$','$p_s$','$f_i$','$m_r$','$m_s$', '$q_r$', '$q_s$', '$t_s$'};

if nargin<3
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
Gampara = [0.878, 6.392, 1, 0.216;%5 fi
    4.645, 2.041, 1, 0.962];%10 ts
% Uniform parameters
% min, max
Unipara =[0.01, 0.03]; %6 mr
end


%Calculate piecewise linear pdf parameters
fa1 = @(l,n) 0.025*2/((l-n)*l);
fa4 = @(m,u) -2*0.025/((m-u)^2);
fb4 = @(m,a4) -m*a4;
fa2 = @(b,u,l,m,a1,a4) 1/((u-l)*(b-l)) * (0.95*2 - a4*(u-b)*(u-m) - a1*l*(b-l)) - a1*l/(b-l);
fb2 = @(a1,a2,l) (a1-a2)*l;
fa3 = @(u,b,a4,b4,a2,b2) 1/(u-b)*(a4*u + b4 - b2 - a2*b);
fb3 = @(a2,b,b2,a3) a2*b + b2 - a3*b;


%% Peicewise Linear
if   ismember(gidx,[1:3,6:8])
    xmin = PLpara(gidx,1);
    xlow = PLpara(gidx,2);
    xbase = PLpara(gidx,3);
    xupp = PLpara(gidx,4);
    xmax = PLpara(gidx,5);

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
    Xpdf = fpiecewise(xgrid);

end

%% Gamma
% X4-fi
if gidx == 4
    a = Gampara(1,1);
    alpha = Gampara(1,2);
    beta = Gampara(1,3);
    b = Gampara(1,4);
    fgamma = @(y) gampdf((y./a).^(1/b),alpha,beta).* (1/(a*b)).*(y./a).^(1/b-1);
    %testpdf=integral(fgamma,0,1.8);
    Xpdf = fgamma(xgrid);

end

% X9, ts
if gidx == 9
    a = Gampara(2,1);
    alpha = Gampara(2,2);
    beta = Gampara(2,3);
    b = Gampara(2,4);

    fgamma = @(y) gampdf((y./a).^(1/b),alpha,beta).* (1/(a*b)).*(y./a).^(1/b-1);
      Xpdf = fgamma(xgrid);

end
%% Uniform

% X5 - mr
if gidx == 5
    Xpdf = 1/(Unipara(2) -Unipara(1))*ones(size(xgrid));
end