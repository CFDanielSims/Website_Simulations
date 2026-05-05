function s = moc_strang_from_M(Mi, M_exit, strang_exit, pos_or_neg_char)

if strcmpi(pos_or_neg_char, 'neg')||strcmpi(pos_or_neg_char, 'negative')||strcmpi(pos_or_neg_char, -1);
    var = 1;
else
    var = -1; 
end

v1 = prandtlmeyer_angled(M_exit);
v2 = prandtlmeyer_angled(Mi);

s = var .* (v1 - v2) + strang_exit;