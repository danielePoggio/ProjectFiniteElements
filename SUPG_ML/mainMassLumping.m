clear all
close all
clc

%% Problema differenziale
u = @(x,y) 2 - 1/(tanh(x + 30*y + 2));
[gradu, d2u] = calculateDerivate(u);
mu = @(x,y) 10.0e-06;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 1.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,y);
gDi = @(x,y) u(x,y);

%% Valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
Pk = 1;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.02;
errorL2ML = zeros(Ktest,1);
errorH1ML = zeros(Ktest,1);
errorLInfML = zeros(Ktest,1);

errorL2P1 = zeros(Ktest,1);
errorH1P1 = zeros(Ktest,1);
errorLInfP1 = zeros(Ktest,1);
Sevec = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/3;
    end
    geom = TriangolatorP1Di(area);
    Area = [geom.support.TInfo.Area].';
    areaTri(l) = max(Area);
    h = sqrt(areaTri(l));
    mk = 1/24;
    Se = sigma(0,0)*2*area/(6*mu(0,0));
    Sevec(l) = Se;
    % valore soluzione esatta nei punti della triangolazione:
    XY = geom.elements.coordinates;
    x = XY(:,1);
    y = XY(:,2);
    Np = length(x);
    sol = zeros(Np,1);
    for k=1:Np
        sol(k) = u(x(k), y(k));
    end
    % calcolo soluzione con Mass Lumping
    uhML = massLumping(geom, mu, sigma, f, gDi, gNe);
    [errorL2ml, errorH1ml] = errorFunctionOld(geom, u, gradu, uhML, Pk);
    errorInfml = norm(sol - uhML, "inf");
    errorL2ML(l) = errorL2ml;
    errorH1ML(l) = errorH1ml;
    errorLInfML(l) = errorInfml;
    % calcolo soluzione con metodo P1
    [uhP1, condA] = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uhP1, Pk);
    errorInf = norm(sol - uhP1, "inf");
    errorL2P1(l) = errorL2;
    errorH1P1(l) = errorH1;
    errorLInfP1(l) = errorInf;
end
figure(1)
loglog(sqrt(areaTri), errorL2ML)
hold on
loglog(sqrt(areaTri), errorL2P1)
legend("Errore Mass Lumping", "Errore P1")
xlabel("h")
ylabel("Errore")
title("Andamento errore norma L2")

figure(2)
loglog(sqrt(areaTri), errorH1ML)
hold on
loglog(sqrt(areaTri), errorH1P1)
legend("Errore Mass Lumping", "Errore P1")
xlabel("h")
ylabel("Errore")
title("Andamento errore norma H1")

figure(3)
loglog(sqrt(areaTri), errorLInfML)
hold on
loglog(sqrt(areaTri), errorLInfP1)
legend("Errore Mass Lumping", "Errore P1")
xlabel("h")
ylabel("Errore")
title("Andamento errore norma L infinito")

pL2ML = polyfit(log(sqrt(areaTri)), log(errorL2ML), 1);

pH1ML = polyfit(log(sqrt(areaTri)), log(errorH1ML), 1);

pL2P1 = polyfit(log(sqrt(areaTri)), log(errorL2P1), 1);

pH1P1 = polyfit(log(sqrt(areaTri)), log(errorH1P1), 1);
