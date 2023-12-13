function w = nodiQuadratura1D(Nq)
%nodiQuadratura1D calcola i pesi per nodi di quadratura equispaziati tra
%[0,1]. Non usare con Nq troppo grande
t = linspace(0,1,Nq);
V = zeros(Nq,Nq);
y = zeros(Nq,1);
for i=0:(Nq-1)
    V(i+1,:) = t.^i;
    y(i+1) = 1/(i+1);
end
w = V\y;
end