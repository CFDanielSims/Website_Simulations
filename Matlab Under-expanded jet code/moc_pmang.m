function s = moc_pmang(v1, phi1, phi2, pos_or_neg_char)

if strcmpi(pos_or_neg_char, 'neg')||strcmpi(pos_or_neg_char, 'negative')||strcmpi(pos_or_neg_char, -1)
    var = 1;
else
    var = -1; 
end

s = var .* (phi1 - phi2) + v1;



%%Returns stream angle on points 2 to the x axis on the same 
%% gamma +/- characteristic when prandtlmeyer angle is known on both positions
%%v1 and v2 and the stream angle is known at the first position phi1

%%Returns in degrees




%%Example data to use in command window

% v2 = prandtlmeyer_angle(2);          %First characteristic
% v1 = prandtlmeyer_angle(2.4436);     %Last characteristic
% phi1 = 11.4137;                      %Angle of starting streamline
% pos_or_neg_char = 'positive';        %Positive or negative characteristic
% 
% moc_strang(v1, v2, phi1, pos_or_neg_char)
