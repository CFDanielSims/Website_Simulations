% expfanslope_from_Mtest

M_exit = 2;
Mi = M_exit:0.01:2.44;
strangle_exit = 0;

s = expfanslope_from_M(2.443, M_exit, strangle_exit, 'positive')

 % view_data_in_table(Mi, s)