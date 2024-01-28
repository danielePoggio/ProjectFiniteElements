clear all
close all
clc

%% Problema differenziale
g = @(x,a) 1./(1 + exp(a*(x-0.5)));
u = @(x,y) g(x^2+y^2,10);
[gradu, d2u] = calculateDerivate(u);
mu = @(x,y) 1.0e-5;
beta = @(x,y) [1,1];
sigma = @(x,y) 0.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,0);
gDi = @(x,y) u(x,y);

%% Valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
Pk = 2;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorL2SUPG = zeros(Ktest,1);
errorH1SUPG = zeros(Ktest,1);
errorL2P2 = zeros(Ktest,1);
errorH1P2 = zeros(Ktest,1);
PeSUPG = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end
    geom = TriangolatorP2Di(area);
    close all
    Area = [geom.support.TInfo.Area].';
    areaTri(l) = max(Area);
    h = sqrt(areaTri(l));
    mk = 1/24;
    Pe = mk*(norm(beta(0,0),2)*h)/(2*mu(0,0));
    PeSUPG(l) = Pe;
    % soluzione con SUPG
    uhSUPG = SUPG(geom, mu, beta, f, gDi, gNe, Pk);
    [errorL2S, errorH1S] = errorFunctionOld(geom, u, gradu, uhSUPG, Pk);
    errorL2SUPG(l) = errorL2S;
    errorH1SUPG(l) = errorH1S;
    % soluzione con P2 classici
    uhP2 = FEMDiNeP2(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uhP2, Pk);
    errorL2P2(l) = errorL2;
    errorH1P2(l) = errorH1;
end
%% calcolo ordini di convergenza

pL2SUPG = polyfit(log(sqrt(areaTri)), log(errorL2SUPG), 1);

pH1SUPG = polyfit(log(sqrt(areaTri)), log(errorH1SUPG), 1);

pL2P2 = polyfit(log(sqrt(areaTri)), log(errorL2P2), 1);

pH1P2 = polyfit(log(sqrt(areaTri)), log(errorH1P2), 1);

figure(1)
loglog(sqrt(areaTri), errorL2P2)
hold on
loglog(sqrt(areaTri), errorL2SUPG)
legend("Errore L2 P2", "Errore L2 SUPG")
xlabel("Numero di Pechlet")
ylabel("Errore in norma L2")