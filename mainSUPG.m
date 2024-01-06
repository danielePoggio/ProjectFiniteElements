clear all
close all
clc

%% Problema differenziale
u = @(x,y) funzione_specifica(x,y);
run("calculateDerivate.m")
gradu = @(x,y) gradu(x,y)';
d2u = @(x,y) [1,0]*Hu(x,y)*[1,0]'+ [0,1]*Hu(x,y)*[0,1]';
mu = @(x,y) 1;
beta = @(x,y) [1.0, 0.0];
sigma = @(x,y) 0.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,0);
gDi = @(x,y) 0;

%% Valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
Pk = 1;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorL2vec = zeros(Ktest,1);
errorH1vec = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end
    geom = TriangolatorP2(area,Pk);
    close all
    Area = [geom.support.TInfo.Area].';
    areaTri(l) = max(Area);
%     uh = SUPG(geom, Pk, mu, beta, f, gDi, gNe);
    uh = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorL2vec(l) = errorL2;
    errorH1vec(l) = errorH1;
end
figure(1)
loglog(sqrt(areaTri), errorL2vec)
title("Andamento errore norma L2")

pL2 = polyfit(log(sqrt(areaTri)), log(errorL2vec), 1);

figure(2)
loglog(sqrt(areaTri), errorH1vec)
title("Andamento errore norma H1")

pH1 = polyfit(log(sqrt(areaTri)), log(errorH1vec), 1);