clear all
close all
clc

%% Problema differenziale
u = @(x,y) 2 - 1/(tanh(x + 30*y + 2));
% run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\calculateDerivate.m")
[gradu, d2u] = calculateDerivate(u);
% mu = @(x,y) 10.0e-06;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
beta = @(x,y) [1.0, 1.0];
sigma = @(x,y) 0.0;
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,y);
gDi = @(x,y) u(x,y);
%% Definiamo triangolazione
Pk = 2;
area = 0.01;
geom = TriangolatorP2Di(area);
Area = [geom.support.TInfo.Area].';
areaTri = max(Area);
h = sqrt(areaTri);

%% Test SUPG
Ktest = 3;
errorL2SUPG = zeros(Ktest,1);
errorH1SUPG = zeros(Ktest,1);
PevSUPG = zeros(Ktest,1);
epsVec = [10.0e-02, 10.0e-04, 10.0e-06];
l = 1;
for eps=epsVec
    mu = @(x,y) eps;
    f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
     mk = 1/24;
    Pe = mk*(norm(beta(0,0),2)*h)/(2*mu(0,0));
    PevSUPG(l) = Pe;
    uh = SUPG(geom, mu, beta, f, gDi, gNe, Pk);
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorL2SUPG(l) = errorL2;
    errorH1SUPG(l) = errorH1;
    l = l + 1;
end

figure(1)
loglog(epsVec, errorL2SUPG)
title("Andamento errore norma L2 SUPG")

pL2SUPG = polyfit(epsVec, log(errorL2SUPG), 1);

figure(2)
loglog(epsVec, errorH1SUPG)
title("Andamento errore norma H1 SUPG")

pH1SUPG = polyfit(epsVec, log(errorH1SUPG), 1);

figure(3)
loglog(epsVec, PevSUPG)
title("Andamento numero di Pe SUPG")


%% Test con P2
Ktest = 3;
errorL2P2 = zeros(Ktest,1);
errorH1P2 = zeros(Ktest,1);
PevP2 = zeros(Ktest,1);
l = 1;
for eps=epsVec
    mu = @(x,y) eps;
    f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
     mk = 1/24;
    Pe = mk*(norm(beta(0,0),2)*h)/(2*mu(0,0));
    PevP2(l) = Pe;
    uh = FEMDiNeP2(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorL2P2(l) = errorL2;
    errorH1P2(l) = errorH1;
    l = l + 1;
end
figure(4)
loglog(epsVec, errorL2P2)
title("Andamento errore norma L2 P2")

pL2P2 = polyfit(epsVec, log(errorL2P2), 1);

figure(5)
loglog(epsVec, errorH1P2)
title("Andamento errore norma H1 P2")

pH1P2 = polyfit(epsVec, log(errorH1P2), 1);

figure(6)
loglog(epsVec, PevP2)
title("Andamento numero di Pe P2")

figure(7)
loglog(PevP2, errorL2P2)
hold on
loglog(PevP2, errorL2SUPG)
legend("Errore L2 P2", "Errore L2 SUPG")
xlabel("Numero di Pechlet")
ylabel("Errore in norma L2")

close all
