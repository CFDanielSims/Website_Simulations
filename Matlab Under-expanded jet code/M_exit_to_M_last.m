function s = M_exit_to_M_last(M_exit, patm_pexit_ratio, gamma)

ratio = (1+((gamma-1)/2)*M_exit.^2).^(-gamma/(gamma-1));

pa_p0 = patm_pexit_ratio.*ratio;

M_last = sqrt((2/(gamma-1))*(((1/pa_p0)^((gamma-1)/gamma))-1));

s = M_last;



