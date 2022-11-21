function Q = QALW(d,q,A,mr,ms,trm,tsm,ts) 


Q = (1-d).*q.*A - (1-d).*mr.*trm - (1-d).*ms.*tsm - (1-d).*ts;
end  