function tTable = tTableforP2plot(ele)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
sizeT = size(ele);
Nele = sizeT(1);
tTable = zeros(4*Nele,3);
for e=0:Nele-1
    tTable(4*e + 1,:) = [ele(e+1,1), ele(e+1,4), ele(e+1,6)];
    tTable(4*e + 2,:) = [ele(e+1,4), ele(e+1,2), ele(e+1,5)];
    tTable(4*e + 3,:) = [ele(e+1,6), ele(e+1,5), ele(e+1,3)];
    tTable(4*e + 4,:) = [ele(e+1,4), ele(e+1,5), ele(e+1,6)];
end
end