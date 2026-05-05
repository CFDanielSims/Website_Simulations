function s = moc_intang(v1plus, v2minus, phi1plus, phi2minus, str_or_pm)

s2 = 0.5.*(v2minus + v1plus) + 0.5.*(phi2minus - phi1plus);
s1 = 0.5.*(v2minus - v1plus) + 0.5.*(phi2minus + phi1plus);

if strcmpi(str_or_pm, 'str')||strcmpi(str_or_pm, 'stream')||strcmpi(str_or_pm, 's');
s = s1;
    else
    s = s2; 
end

