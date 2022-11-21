function [y]=PSC_87(x)
% Plante et al. (1987) DA of Treatment options in Pyriform Sinus Carcinoma

ds = x(:,1);
pr = x(:,2);
ps = x(:,3);
fi = x(:,4);
mr = x(:,5);
ms = x(:,6);
qr = x(:,7);
qs = x(:,8);
ts = x(:,9);

Nc = size(x,1);
% dr, Tsm and Trm are fixed at base values
dr = 0.005*ones(Nc,1);%%
tsm = 2.5*ones(Nc,1);% tsc
trm = 2*ones(Nc,1);%% tr
tr = 0*ones(Nc,1);%% 

% alternative 1 radiation therapy
A = Autility(pr,1);
Q1 = QALW(dr,qr,A,mr,0,trm,0,0) ;

% alternative 2 Surgery
A = Autility(ps,2);
Q2 = QALW(ds,qs,A,0,ms,0,tsm,ts) ;

% alternative 3 Pre-operative irradiaition then surgery
A = Autility(fi.*ps,3);
Q3 = QALW(dr+ds -dr.*ds,qr.*qs,A,mr,ms,tr,tsm,ts) ; 

% alternative 4 surgery then post-operative irradiation
A = Autility(fi.*ps,4);
Q4 = QALW(dr+ds -dr.*ds,qr.*qs,A,mr,ms,tr,tsm,ts) ; 

y = [Q1,Q2,Q3,Q4];
end