function S = getSchnittpunkt(curve)
l = length(curve);
off = 5;

initLin = polyfit(1:10,curve(1:10),1);
a1 = initLin(1);
b1 = initLin(2);

endLin = polyfit(l-off:-1:l-off-9, curve(l-off:-1:l-off-9),1);
a2 = endLin(1);
b2 = endLin(2);

Sx = floor((b2 - b1)/(a1 - a2));
S = curve(Sx);

end