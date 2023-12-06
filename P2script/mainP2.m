clear all
close all
clc

%% Eseguo Triangolazione sul Dominio
area = 0.01;
Pk = 2;
geom = TriangolatorP2(area);
close all

%% Problema differenziale
u = @(x,y) 16*x*(1-x)*y*(1-y);
gradu = @(x,y) 16*[y*(1-y)*(1 -2*x), x*(1-x)*(1 - 2*y)]';
mu = @(x,y) 1.0;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(x,y) 32*(x*(1-x) + y*(1-y)); % + 16*(1-2*x)*(y*(1-y));
gNe = @(x,y) -16*x*(1-x);
gDi = @(x,y) 0;

%% Soluzione problema discretizzato
% uh = FEMDiNeP2(geom, mu, beta, sigma, f, gDi, gNe);
uh = SUPG(geom, Pk, mu, beta, f, gDi, gNe);

%% Plot soluzione approssimata
% tTable = tTableforP2plot(geom.elements.triangles);
XY = geom.elements.coordinates;
x = XY(:,1);
y = XY(:,2);
figure(1)
tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
trisurf(tTable, x, y, uh);
title("Grafico funzione approssimata")

% Calcolo errore
% Pk = 2;
% [errorL2, errorH1] = errorFunction(geom, u, gradu, uh, Pk);

%% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
Pk = 2;
Ktest = 4;
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
%     if area == 0
    geom = TriangolatorP2(area);
    close all
    Area = [geom.support.TInfo.Area].';
    maxArea = max(Area);
    areaTri(l) = maxArea;
    uh = FEMDiNeP2(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunction(geom, u, gradu, uh, Pk);
    clear geom
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
