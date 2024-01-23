function j_ = localIndexP3(j)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if any([1,2,3,4] == j)
    j_ = j;
elseif any([5,6] == j)
    j_ = 1;
elseif any([7,8] == j)
    j_ = 2;
elseif any([9, 10] == j)
    j_ = 3;
end
end

