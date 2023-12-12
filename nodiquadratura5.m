q = 5;
t = linspace(0,1,q);
V = zeros(q,q);
b = zeros(q,1);
for i=0:(q-1)
    V(i+1,:) = t.^i;
    b(i+1) = 1/(i+1);
end
w = V\b;

phi1 = @(t) 2*(t-0.5)*(t-1);
phi2 = @(t) -4*t*(t-1);
phi3 = @(t) 2*t*(t-0.5);

phi = @(t) [phi1(t), phi2(t), phi3(t)]';
phiM = @(t) phi(t)*phi(t)';
phiTensor = zeros(3,3,q);
for k=1:q
    phiTensor(:,:,k) = w(k)*phiM(t(k)); 
end
phiMatrix = sum(phiTensor,3);
gN = [0.4 0.7 10]';
x = phiMatrix*gN;


