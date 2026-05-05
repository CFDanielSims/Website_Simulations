function [ycharlines, yjetb] = solver_loop1(Mi, Strangd, Strang_exit, pos_or_neg_char, x, init_char, ax_of_sym, gamma, v_last, radius)

N = length(Mi);
Nv = 7;                                            %Number of variables along a characteristic line 

s = zeros(N, N, Nv);               %N is amount of lines

Mi;                                              %Mach angles from 2 to 2.443
Strangd;                                            %Initial streamline angle (Nx1) (Strangd)
% phi2;                                            %Axis of symmetry (Always 0)    (Strang_exit)
v1 = prandtlmeyer_angled(Mi);                                %Initial Prandtl Meyer angle (Nx1)

Nx = length(x);
N = length(Mi);

init_char; %Initial characteristic lines, formed in master code, Nx by N

char_lines = zeros(N, Nx);
x_int = zeros(N);
y_int = zeros(N);
x_index = zeros(N);

%1st index is Prandtl Meyer angle, 
%2nd index is Stream Angle
%3rd index is x coordinate intersection
%4th index is y coordinate intersection
%5th index is gamma plus slope characteristic
%6th index is gamma minus slope characteristic
%7th index is index of x at which x_int coordinate is found closest to x

for j = 1:N
    
    u = moc_pmang_from_M(Mi(j), Strangd(j), Strang_exit, pos_or_neg_char);  %Store Prandtl Meyer angle in each of the diagonals
    s(j, j, 1) = u;

        for i = j+1:N
            l = moc_intang(s(j, j, 1), v1(i), s(j, j, 2), Strangd(i), 'pm');
            k = moc_intang(s(j, j, 1), v1(i), s(j, j, 2), Strangd(i), 'str');

            s(i, j, 1) = l;                                      % Store the result of the first integration
            s(i, j, 2) = k;                                      % Store the result of the second integration
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

    [x_int(j, j), y_int(j, j), x_index(j, j)] = intersection_finder_index(x, x, char_lines(j, :), ax_of_sym); %Give intersection points with axis of symmetry

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

% view_data_in_table_array(Strangd);




%JET BOUNDARY LINE%%%%%%%%%%%%%%%%%

jetb = zeros(N+1, 4);               %Jet boundary information for each point

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

% view_data_in_table_array(ycharlines);

yjetb = zeros(1, Nx);

j = 1;
[~, f] = line_through_point(jetb(j, 2), jetb(j, 3), tand(jetb(j, 1)));

yjetb = f(x);

% j = 3;
% view_data_in_table_array(yjetb(s(N, j, 7):end));

x_intjb = zeros(1, N);
y_intjb = zeros(1, N);
x_indexjb = zeros(1, N);

for j = 2:N+1
    
    % view_data_in_table_array(s(N, j, 7))
    % view_data_in_table(yjetb(s(N, j, 7):end), char_lines(j, s(N, j, 7):end));
    [x_intjb(j), y_intjb(j), x_indexjb(j)] = intersection_finder_index(x, x(s(N, j-1, 7):end), yjetb(s(N, j-1, 7):end), char_lines(j-1, s(N, j-1, 7):end));
    
    jetb(j, 2) = x_intjb(j);
    
    t1 = x_int(j);
    t2 = x(s(N, j-1, 7):end);
    t3 = yjetb(s(N, j-1, 7):end);
    t4 = char_lines(j-1, s(N, j-1, 7):end);

    jetb(j, 3) = y_intjb(j);
    jetb(j, 4) = x_indexjb(j);

    [~, f] = line_through_point(jetb(j, 2), jetb(j, 3), tand(jetb(j, 1)));
    
    yjetb(jetb(j, 4):end) = f(x(jetb(j, 4):end));

end
% view_data_in_table_array(yjetb)


%BOUNCE FROM JET BOUNDARY

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

% viewarray(s2(:, :, 1), s2(:, :, 2));
for j = 1:N
    for i = 1:N
    s2(i, j, 5) = expfanslope_from_pmstrang(s2(i, j, 1), s2(i, j, 2), 'pos', gamma);   %Positive slopes
    s2(i, j, 6) = expfanslope_from_pmstrang(s2(i, j, 1), s2(i, j, 2), 'neg', gamma);   %Negative slopes
    end
end

% j = 5;
% plot(x(s(N, j, 7):end), char_lines(j, s(N, j, 7):end), 'r', x(s(N, j, 7):end), yjetb(s(N, j, 7):end), 'b');

for j = 1:N

    [x_int2(j, j), y_int2(j, j), x_index2(j, j)] = intersection_finder_index(x, x(s(N, j, 7):end), char_lines(j, s(N, j, 7):end), yjetb(s(N, j, 7):end)); %Give intersection points with jet boundary

    s2(j, j, 3) = x_int2(j, j);      %x coordinate
    s2(j, j, 4) = y_int2(j, j);      %y coordinate
    s2(j, j, 7) = x_index2(j, j);    %x index

    [~, f] = line_through_point(s2(j, j, 3), s2(j, j, 4), s2(j, j, 6));

    char_lines(j, s2(j, j, 7):end) = f(x(s2(j, j, 7):end));
    % view_data_in_table_array(s2(:, :, 7));
    

    for i = j+1:N
        
        % view_data_in_table_array(x(s2(i-1, j, 7):end))
        % view_data_in_table_array(s2(i-1, j, 7))
        % plot(x(s2(i-1, j, 7):end), char_lines(j, s2(i-1, j, 7):end), 'b', x(s2(i-1, j, 7):end), char_lines(i, s2(i-1, j, 7):end), 'r')
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

viewarray(s(:, :, 1), s(:, :, 2));


%2ND BOUNCE FROM SYMMETRY LINE

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

% viewarray(s2(:, :, 1), s2(:, :, 2), s(:, :, 1), s(:, :, 2));

for j = 1:N
    s3(j, j, 2) = 0;
    % viewarray(s2(N, j, 1), s2(N, j, 2), s3(j, j, 2));
    u = moc_pmang(s2(N, j, 1), s2(N, j, 2), s3(j, j, 2), pos_or_neg_char);  %Store Prandtl Meyer angle in each of the diagonals
    s3(j, j, 1) = u;
    
        for i = j+1:N
            
            l = moc_intang(s3(j, j, 1), s2(N, i, 1), s3(j, j, 2), s2(N, i, 2), 'pm'); %First pm and strang has to be from the gamma + characteristic
            k = moc_intang(s3(j, j, 1), s2(N, i, 1), s3(j, j, 2), s2(N, i, 2), 'str');

            t1 = s3(j, j, 1);
            t2 = s2(N, i, 1);
            t3 = s3(j, j, 2);
            t4 = s2(N, i, 2);

            s3(i, j, 1) = l;                                      % Store the result of the first integration
            s3(i, j, 2) = k;                                      % Store the result of the second integration
        end
end

% viewarray(s3(:, :, 1), s3(:, :, 2));
% viewarray(s(:, :, 1), s(:, :, 2));
% viewarray(v1(:), s2(N, :, 1), Strangd(:), s2(N, :, 2));


for j = 1:N
    for i = 1:N
    s3(i, j, 5) = expfanslope_from_pmstrang(s3(i, j, 1), s3(i, j, 2), 'pos', gamma);   %Positive slopes
    s3(i, j, 6) = expfanslope_from_pmstrang(s3(i, j, 1), s3(i, j, 2), 'neg', gamma);   %Negative slopes
    end
end

% viewarray(s3(:, :, 5), s3(:, :, 6))
% j = 1;
% plot(x(s2(N, j, 7):end), char_lines(j, s2(N, j, 7):end), 'r', x(s2(N, j, 7):end), yjetb(s2(N, j, 7):end), 'b');

for j = 1:N

    [x_int3(j, j), y_int3(j, j), x_index3(j, j)] = intersection_finder_index(x, x(s2(N, j, 7):end), char_lines(j, s2(N, j, 7):end), ax_of_sym(s2(N, j, 7):end)); %Give intersection points with jet boundary

    s3(j, j, 3) = x_int3(j, j);      %x coordinate
    s3(j, j, 4) = y_int3(j, j);      %y coordinate
    s3(j, j, 7) = x_index3(j, j);    %x index

    [~, f] = line_through_point(s3(j, j, 3), s3(j, j, 4), s3(j, j, 5));

    char_lines(j, s3(j, j, 7):end) = f(x(s3(j, j, 7):end));
    % view_data_in_table_array(s3(:, :, 7));
    for i = j+1:N
        
        % view_data_in_table_array(x(s3(i-1, j, 7):end))
        % view_data_in_table_array(s3(i-1, j, 7))
        % plot(x(s3(i-1, j, 7):end), char_lines(j, s3(i-1, j, 7):end), 'b', x(s3(i-1, j, 7):end), char_lines(i, s3(i-1, j, 7):end), 'r')
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

% viewarray(s(:, :, 1), s2(:, :, 1), s3(:, :, 1));









% viewarray(s2(:, :, 1));
% viewarray(s2(:, :, 2));

ycharlines = char_lines;


















% view_data_in_table_array(jetb); %%Why is index not increasing???
% plot(x, yjetb, 'w')



    



% view_data_in_table_array(jetb);



% 
% t = s(j, j, 7);
%     t25 = ybef:s(j, j, 7);
%     t2 = x(ybef:s(j, j, 7));
%     t3 = char_lines(j, ybef:s(j, j, 7));
%     t4 = f(x(ybef:s(j, j, 7)));

 % if j>1
 %        ybef = s(j-1, j-1, 7);
 %    else
 %        ybef = 1;
 %    end

     % view_data_in_table_array(char_lines);

       % [x_int(j, j), y_int(j, j)] = intersection_finder_index(x, char_lines(:, i), char_lines()