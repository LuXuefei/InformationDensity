function [infogain,astarstar,X,Y]  = EVPI2ddouleloopgrid(xgrid,ygrid,gidx,astar,Nc,PLpara,Gampara,Unipara)


if nargin<5
    Nc = 10000;
   %Piecewise linear parameters
% Min, Lower, Base, Upper, Max

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

[X,Y] = meshgrid(xgrid,ygrid);


M= length(xgrid);
infogain=zeros(M,M);
astarstar = zeros(M,M);
EYAstarstarXu = zeros(M,M); % EUastarstar_Xu
EYAstarXu = zeros(M,M); % EUastar_Xu

x = InputGenerator_Cancer(Nc,PLpara,Gampara,Unipara);

%for m=1:M
for mi = 1:M % column
    for mj = 1:M

        % for k=1:length(xgrid)
        xu = x;
        xu(:,gidx(1)) = X(mi,mj);
        xu(:,gidx(2)) = Y(mi,mj);

        [yaxu]=PSC_87(xu);

%                % Hazen's way -- the same as below
%        [EYAstarstarXu(mi,mj),astarstar(mi,mj)]  =  max(mean(yaxu),[],2);
%        EYAstarXu(mi,mj) = mean(yaxu(:,astar));

        EYaXu = mean(yaxu);
        [EYAstarstarXu(mi,mj),astarstar(mi,mj)] = max(EYaXu);
        EYAstarXu(mi,mj) = EYaXu(astar);
        % end
        infogain(mi,mj) = EYAstarstarXu(mi,mj) - EYAstarXu(mi,mj);

    end
end

% EVPI= mean(EYAstarstarXu - EYAstarXu);%not correct
end

