function s = expfanslope_from_pmstrang(pm, strang, pos_or_neg_char, gamma)

 if strcmpi(pos_or_neg_char, 'neg')||strcmpi(pos_or_neg_char, 'negative')||strcmpi(pos_or_neg_char, -1)
        var = -1;
    else
        var = 1; 
 end

M = M_from_pm(pm, gamma);
s = tand(strang + var .* asind(1/M));
