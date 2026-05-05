function s = moc_pmang_from_M(Mi, Strangd, Strang_exit, pos_or_neg_char);
N = length(Mi);

phi1 = Strangd;
phi2 = Strang_exit .* ones(N, 1);

if strcmpi(pos_or_neg_char, 'neg')||strcmpi(pos_or_neg_char, 'negative')||strcmpi(pos_or_neg_char, -1);
    var = 1;
else
    var = -1; 
end

v1 = prandtlmeyer_angled(Mi);

s = var .* (phi1 - phi2) + v1;




%%Returns Prandtl Meyer angle on points 2 on the same 
%% gamma +/- characteristic when stream angle is known on both positions
%%phi1 and phi2 and the Prandtl Meyer angle is known at the first position
%%v1

%%Returns in degrees

%phi1 is Strangd
%phi2 is Strang_exit = 0
%v1 is initial prandtlmeyer angle





%%Example data to use in command window

% v2 = prandtlmeyer_angled(2);          %First characteristic
% v1 = prandtlmeyer_angled(2.4436);     %Last characteristic
% phi1 = 11.4137;                      %Angle of starting streamline
% pos_or_neg_char = 'positive';        %Positive or negative characteristic
% 
% moc_strang(v1, v2, phi1, pos_or_neg_char)
