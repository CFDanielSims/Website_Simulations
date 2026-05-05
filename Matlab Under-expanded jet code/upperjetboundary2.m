function u = upperjetboundary2(x_int, y_int, x, phi, pm1, pp0ratio, n1)
    
    M2 = M_from_statstag_ratio(pp0ratio);
    
    pm2 = prandtlmeyer_angled(M2);

    j = phi + tand(pm1 - pm2)
   
    [b, func] = line_through_point(x_int, y_int, j);

    jet1b1 = func(x);

    u = jet1b1(n1:end);