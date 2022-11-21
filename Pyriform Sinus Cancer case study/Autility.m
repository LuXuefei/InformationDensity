function [ua] = Autility (pii,type)


% Table 1 Plante 1987

if type == 1 % primary radiation
ua = pii.*( 485) + (1-pii).* ( 0.8.*( 24 )  + 0.15.*(72) + 0.05.* (120)  );
elseif type  == 2 % primary surgery
  ua = pii.*( 485) + (1-pii).* ( 0.7.*( 24.)  + 0.18.*(72) + 0.12.* (120 )  );
elseif type == 3 % preoperative irradiation, then surgery
ua = pii.*( 485) + (1-pii).* ( 0.73.*( 24)  + 0.25.*(72) + 0.02.* (120)  );
elseif type == 4 % surgery then post-operative irradiation
ua = pii.*( 485 ) + (1-pii).* ( 0.70.*( 24)  + 0.06.*(72) + 0.24.* (120 )  );
end

end
