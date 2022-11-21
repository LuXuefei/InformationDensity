function EVPI = EVPIdouleloop2dim(x,indx)

[n,ni]  = size(x);

% unconditional
[ya]=PSC_87(x);
% a^*:
EYa = mean(ya,1);
[MEYa,astar] = max(EYa);

% this is to replicate Table 1,Felli Hazen 1999 HE
%xb = repmat([ 0.05, 0.1,0.25, 1.3,  0.01,0.3,0.9,0.7, 5],n,1);% paper baseline
%xb = repmat([0.0575 0.0797 0.10 0.3518 1.2952 0.02  0.3089  0.8627 0.7302 9.1770],n,1); % real mean\
%real mean Unifrom - result diff from Table6.!!



%for i=1:ni %input
    xgrid = x(:,indx);
    astarstar = nan(n,1);
    EYAstarstarXu = nan(n,1);
    EYAstarXu = nan(n,1);


    for k=1:length(xgrid)
        xu = x;
        %xu = xb; % testbaseline
        xu(:,indx) = ones(n,1)*xgrid(k,:);
        [yaxu]=PSC_87(xu);

       % Hazen's way -- the same as below
       EYAstarstarXu(k) =  max(mean(yaxu),[],2);
       EYAstarXu(k) = mean(yaxu(:,astar));

%        EYaXu = mean(yaxu);
%        [EYAstarstarXu(k),astarstar(k)] = max(EYaXu);
%        EYAstarXu(k) = EYaXu(astar);
    end
    infogain = EYAstarstarXu - EYAstarXu;
    EVPI = mean(infogain);

%     [Xpdf1] = PSC_87_density(indx(1),x(:,indx(1)));
%     [Xpdf2] = PSC_87_density(indx(2),x(:,indx(2)));
% 
%     minx = min(x);
% maxx = max(x);
% V = (maxx(indx(1)) - minx(indx(1)) ) * (maxx(indx(2)) - minx(indx(2)) );
% 
%    
%     
%     EVPI = V*mean( Xpdf1 .* Xpdf2 .* infogain);



end
%end

