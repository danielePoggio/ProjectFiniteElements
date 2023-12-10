function w = calcoloNeP2nodi()
A = [[1,1,1];[0 0.5 1];[0 0.25 1]];
b = [1 0.5 1/3]';
w = A\b;
end