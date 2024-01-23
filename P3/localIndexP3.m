function j_ = localIndexP3(j)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if any([1,4,7,10] == j)
    j_ = 1;
elseif any([2,5,8] == j)
    j_ = 2;
elseif any([3,6,9] == j)
    j_ = 3;
end
end

