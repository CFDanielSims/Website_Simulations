function [ycharlines, yjetb, s, s2, s3] = solver_loop1clean(Mi, Strangd, Strang_exit, pos_or_neg_char, x, init_char, ax_of_sym, gamma, v_last, radius)

N = length(Mi);
Nv = 7;                                            %Number of variables along a characteristic line 

s = zeros(N, N, Nv);                                %N is amount of lines
v1 = prandtlmeyer_angled(Mi);                       %Initial Prandtl Meyer angle (Nx1)

Nx = length(x);                                     %Number of x coordinates
N = length(Mi);                                     %Number of characteristic lines

x_int = zeros(1, N);
y_int = zeros(1, N);
x_index = zeros(1, N);

%FIRST BOUNCING LINES FROM SYMMETRY LINE

%1st index is Prandtl Meyer angle, 
%2nd index is Stream Angle
%3rd index is x coordinate intersection
%4th index is y coordinate intersection
%5th index is gamma plus slope characteristic
%6th index is gamma minus slope characteristic
%7th index is index of x at which x_int coordinate is found closest to x

for j = 1:N
    
    %Store Prandtl Meyer angle in each of the diagonals
    u = moc_pmang_from_M(Mi(j), Strangd(j), Strang_exit, pos_or_neg_char);  
    s(j, j, 1) = u;

        for i = j+1:N
            l = moc_intang(s(j, j, 1), v1(i), s(j, j, 2), Strangd(i), 'pm');
            k = moc_intang(s(j, j, 1), v1(i), s(j, j, 2), Strangd(i), 'str');

            s(i, j, 1) = l;    
            s(i, j, 2) = k;    
        end
end

for j = 1:N
    for i = 1:N
    s(i, j, 5) = expfanslope_from_pmstrang(s(i, j, 1), s(i, j, 2), 'pos', gamma);   %Positive slopes
    s(i, j, 6) = expfanslope_from_pmstrang(s(i, j, 1), s(i, j, 2), 'neg', gamma);   %Negative slopes
    end
end

for j = 1:N
    if j == 1
        char_lines = init_char;
    end
    
    %Give intersection points with axis of symmetry
    [x_int(j, j), y_int(j, j), x_index(j, j)] = intersection_finder_index(x, x, char_lines(j, :), ax_of_sym); 

    s(j, j, 3) = x_int(j, j);
    s(j, j, 4) = y_int(j, j);
    s(j, j, 7) = x_index(j, j);

    [~, f] = line_through_point(s(j, j, 3), s(j, j, 4), s(j, j, 5));

    char_lines(j, s(j, j, 7):end) = f(x(s(j, j, 7):end));

    for i = j+1:N
        [x_int(i, j), y_int(i, j), x_index(i, j)] = intersection_finder_index(x, x(s(i-1, j, 7):end), char_lines(j, s(i-1, j, 7):end), char_lines(i, s(i-1, j, 7):end)); %Finds the intersection point between the l1, 1 and l2

        s(i, j, 3) = x_int(i, j);
        s(i, j, 4) = y_int(i, j);
        s(i, j, 7) = x_index(i, j);

        [~, fplus] = line_through_point(s(i, j, 3), s(i, j, 4), s(i, j, 5));
        char_lines(j, s(i, j, 7):end) = fplus(x(s(i, j, 7):end));                           %Line from second intersection point to end (First intersection point being intersection with y = 0)

        [~, fmin] = line_through_point(s(i, j, 3), s(i, j, 4), s(i, j, 6));
        char_lines(i, s(i, j, 7):end) = fmin(x(s(i, j, 7):end));
    end
end

%JET BOUNDARY LINE
jetb = zeros(N+1, 4);               %Initial jet boundary information for each point

%1st column is angle
%2nd column is x_intersection coordinate
%3rd column is y_intersection coordinate
%4th column is x_column of intersection

for j = 1:N+1

    if j == 1
        jetb(1, 1) = Strangd(N);
    else
        jetb(j, 1) = moc_strang(s(N, j-1, 1), v_last, s(N, j-1, 2), 'pos');
    end
end

jetb(1, 3) = radius;                 %Starting point for the jet boundary

j = 1;
[~, f] = line_through_point(jetb(j, 2), jetb(j, 3), tand(jetb(j, 1)));

yjetb = f(x);

x_intjb = zeros(1, N);
y_intjb = zeros(1, N);
x_indexjb = zeros(1, N);

for j = 2:N+1
    
    [x_intjb(j), y_intjb(j), x_indexjb(j)] = intersection_finder_index(x, x(s(N, j-1, 7):end), yjetb(s(N, j-1, 7):end), char_lines(j-1, s(N, j-1, 7):end));
    
    jetb(j, 2) = x_intjb(j);
    jetb(j, 3) = y_intjb(j);
    jetb(j, 4) = x_indexjb(j);

    [~, f] = line_through_point(jetb(j, 2), jetb(j, 3), tand(jetb(j, 1)));
    
    yjetb(jetb(j, 4):end) = f(x(jetb(j, 4):end));

end

%BOUNCING LINES FROM JET BOUNDARY

%1st index is Prandtl Meyer angle, 
%2nd index is Stream Angle
%3rd index is x coordinate intersection
%4th index is y coordinate intersection
%5th index is gamma plus slope characteristic
%6th index is gamma minus slope characteristic
%7th index is index of x at which x_int coordinate is found closest to x

x_int2 = zeros(N);
y_int2 = zeros(N);
x_index2 = zeros(N);
s2 = zeros(N, N, Nv);

for j = 1:N

    u = v_last;  %Store Prandtl Meyer angle in each of the diagonals, stays consistent over the jet boundary
    p = jetb(j+1, 1); %Strang of jet boundary in pointinfo matrix s2
    s2(j, j, 1) = u;
    s2(j, j, 2) = p;

        for i = j+1:N
            
            l = moc_intang(s(N, i, 1), s2(j, j, 1), s(N, i, 2), s2(j, j, 2), 'pm'); %First pm and strang has to be from the gamma + characteristic
            k = moc_intang(s(N, i, 1), s2(j, j, 1), s(N, i, 2), s2(j, j, 2), 'str');

            s2(i, j, 1) = l;                                      % Store the result of the first integration
            s2(i, j, 2) = k;                                      % Store the result of the second integration
        end
end

for j = 1:N
    for i = 1:N
    s2(i, j, 5) = expfanslope_from_pmstrang(s2(i, j, 1), s2(i, j, 2), 'pos', gamma);   %Positive slopes
    s2(i, j, 6) = expfanslope_from_pmstrang(s2(i, j, 1), s2(i, j, 2), 'neg', gamma);   %Negative slopes
    end
end

for j = 1:N

    [x_int2(j, j), y_int2(j, j), x_index2(j, j)] = intersection_finder_index(x, x(s(N, j, 7):end), char_lines(j, s(N, j, 7):end), yjetb(s(N, j, 7):end)); %Give intersection points with jet boundary

    s2(j, j, 3) = x_int2(j, j);      %x coordinate
    s2(j, j, 4) = y_int2(j, j);      %y coordinate
    s2(j, j, 7) = x_index2(j, j);    %x index

    [~, f] = line_through_point(s2(j, j, 3), s2(j, j, 4), s2(j, j, 6));

    char_lines(j, s2(j, j, 7):end) = f(x(s2(j, j, 7):end));

    for i = j+1:N
       
        [x_int2(i, j), y_int2(i, j), x_index2(i, j)] = intersection_finder_index(x, x(s2(i-1, j, 7):end), char_lines(j, s2(i-1, j, 7):end), char_lines(i, s2(i-1, j, 7):end)); %Finds the intersection point between the l1, 1 and l2

        s2(i, j, 3) = x_int2(i, j);
        s2(i, j, 4) = y_int2(i, j);
        s2(i, j, 7) = x_index2(i, j);

        [~, fplus] = line_through_point(s2(i, j, 3), s2(i, j, 4), s2(i, j, 5));
        char_lines(i, s2(i, j, 7):end) = fplus(x(s2(i, j, 7):end));                           %Line from second intersection point to end (First intersection point being intersection with y = 0)
            
        [~, fmin] = line_through_point(s2(i, j, 3), s2(i, j, 4), s2(i, j, 6));
        char_lines(j, s2(i, j, 7):end) = fmin(x(s2(i, j, 7):end));

    end
end


%2ND BOUNCING LINES FROM SYMMETRY LINE

%1st index is Prandtl Meyer angle, 
%2nd index is Stream Angle
%3rd index is x coordinate intersection
%4th index is y coordinate intersection
%5th index is gamma plus slope characteristic
%6th index is gamma minus slope characteristic
%7th index is index of x at which x_int coordinate is found closest to x

x_int3 = zeros(N);
y_int3 = zeros(N);
x_index3 = zeros(N);
s3 = zeros(N, N, Nv);

for j = 1:N
    s3(j, j, 2) = 0;
    u = moc_pmang(s2(N, j, 1), s2(N, j, 2), s3(j, j, 2), pos_or_neg_char);  %Store Prandtl Meyer angle in each of the diagonals
    s3(j, j, 1) = u;
    
        for i = j+1:N
            
            l = moc_intang(s3(j, j, 1), s2(N, i, 1), s3(j, j, 2), s2(N, i, 2), 'pm'); %First pm and strang has to be from the gamma + characteristic
            k = moc_intang(s3(j, j, 1), s2(N, i, 1), s3(j, j, 2), s2(N, i, 2), 'str');

            s3(i, j, 1) = l;                                      % Store the result of the first integration
            s3(i, j, 2) = k;                                      % Store the result of the second integration
        end
end

for j = 1:N
    for i = 1:N
    s3(i, j, 5) = expfanslope_from_pmstrang(s3(i, j, 1), s3(i, j, 2), 'pos', gamma);   %Positive slopes
    s3(i, j, 6) = expfanslope_from_pmstrang(s3(i, j, 1), s3(i, j, 2), 'neg', gamma);   %Negative slopes
    end
end

for j = 1:N

    [x_int3(j, j), y_int3(j, j), x_index3(j, j)] = intersection_finder_index(x, x(s2(N, j, 7):end), char_lines(j, s2(N, j, 7):end), ax_of_sym(s2(N, j, 7):end)); %Give intersection points with jet boundary

    s3(j, j, 3) = x_int3(j, j);      %x coordinate
    s3(j, j, 4) = y_int3(j, j);      %y coordinate
    s3(j, j, 7) = x_index3(j, j);    %x index

    [~, f] = line_through_point(s3(j, j, 3), s3(j, j, 4), s3(j, j, 5));

    char_lines(j, s3(j, j, 7):end) = f(x(s3(j, j, 7):end));
    
    for i = j+1:N
     
        [x_int3(i, j), y_int3(i, j), x_index3(i, j)] = intersection_finder_index(x, x(s3(i-1, j, 7):end), char_lines(j, s3(i-1, j, 7):end), char_lines(i, s3(i-1, j, 7):end)); %Finds the intersection point between the l1, 1 and l2

        s3(i, j, 3) = x_int3(i, j);
        s3(i, j, 4) = y_int3(i, j);
        s3(i, j, 7) = x_index3(i, j);

        [~, fplus] = line_through_point(s3(i, j, 3), s3(i, j, 4), s3(i, j, 5));
        char_lines(j, s3(i, j, 7):end) = fplus(x(s3(i, j, 7):end));                           %Line from second intersection point to end (First intersection point being intersection with y = 0)
            
        [~, fmin] = line_through_point(s3(i, j, 3), s3(i, j, 4), s3(i, j, 6));
        char_lines(i, s3(i, j, 7):end) = fmin(x(s3(i, j, 7):end));

    end
end

ycharlines = char_lines;



















    


