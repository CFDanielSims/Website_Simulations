v2 = prandtlmeyer_angled(2.4436);          %First characteristic
v1 = prandtlmeyer_angled(2.4436);     %Last characteristic
phi1 = -11.4137;                     %Angle of starting streamline
pos_or_neg_char = 'positive';        %Positive or negative characteristic

mu = moc_strang(v1, v2, phi1, pos_or_neg_char)