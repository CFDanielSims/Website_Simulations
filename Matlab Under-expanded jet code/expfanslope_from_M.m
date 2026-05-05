function s = expfanslope_from_M(Mi, M_exit, strangle_exit, pos_or_neg_char)
    if strcmpi(pos_or_neg_char, 'neg')||strcmpi(pos_or_neg_char, 'negative')||strcmpi(pos_or_neg_char, -1);
        var = 1;
    else
        var = -1; 
    end


    s = tand(var.*(-prandtlmeyer_angled(Mi) + prandtlmeyer_angled(M_exit)) + strangle_exit - asind(1./Mi));

     % view_data_in_table(Mi, asind(1./Mi), -prandtlmeyer_angled(Mi), prandtlmeyer_angled(M_exit), s)

%Only used to determine the angle of the slopes from the initial
%characteristic lines