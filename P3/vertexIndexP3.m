function j_ = vertexIndexP3(j)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if any([1,2,3] == j)
    j_ = 1;
elseif any([4,5,6] == j)
    j_ = 2;
elseif any([7,8,9] == j)
    j_ = 3;
elseif 10 == j
    j_ = 4;
end
end

