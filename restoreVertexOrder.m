function [outputArg1,outputArg2] = restoreVertexOrder(ele, XY)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sizeT = size(ele);
Nele = sizeT(1);
Nv = sizeT(2);
for e=1:Nele
    idxVertex = ele(e,:);
    % vertici del triangolo
    p1 = XY(ele(e,1),:);
    p2 = XY(ele(e,2),:);
    p3 = XY(ele(e,3),:);
    
    for idx = 4:Nv

    end

end
end