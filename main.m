clear all
close all
clc

%% Eseguo Triangolazione sul Dominio
area = 0.02;
geom = Triangolator(area);
close all

% %% Assemblaggio FEM + risoluzione sistema lineare
% mu = @(x,y) 1.0;
% beta = @(x,y) [0.0, 0.0];
% sigma = @(x,y) 0.0;
% f = @(x,y) 32*(x*(1-x) + y*(1-y));
% gNe = @(x,y) -16*x*(1-x);
% gDi = @(x,y) 0;

%% Altro problema
mu = @(x,y) 1.0;
beta = @(x,y) [1.0, 0.0];
sigma = @(x,y) 1.0;
f = @(x,y) 8*y*(1-y)*(2-3*x-x^2) + 8*x*(1-x)*(4+y-y^2);
gNe = @(x,y) -16*x*(1-x);
gDi = @(x,y) 0;
uh = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
%% Plot soluzione approssimata
XY = geom.elements.coordinates;
x = XY(:,1);
y = XY(:,2);
figure(1)
tri = delaunay(x, y); % Genera la matrice di connettività dei triangoli
trisurf(tri, x, y, uh);
title("Grafico funzione approssimata")

%% calcoliamo stima dell'errore
u = @(x,y) 16*x*(1-x)*y*(1-y);
gradu = @(x,y) 16*[y*(1-y)*(1 -2*x), x*(1-x)*(1 - 2*y)]';
[errorL2, errorH1] = errorFunction(geom, u, gradu, uh);
% [errorL2quad, errorH1quad] = errorFunction(geom, u, gradu, uh_quad);
%% Andiamo a vedere come estrarre il valore minimo e massimo dell'area nella triangolazione
Area = [geom.support.TInfo.Area].';
areaMax = max(Area);

%% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
mu = @(x,y) 1.0;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(x,y) 32*(x*(1-x) + y*(1-y));
gNe = @(x,y) -16*x*(1-x);
gDi = @(x,y) 0;
u = @(x,y) 16*x*(1-x)*y*(1-y);
gradu = @(x,y) 16*[y*(1-y)*(1 -2*x), x*(1-x)*(1 - 2*y)]';
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
    geom = Triangolator(area);
    close all
    Area = [geom.support.TInfo.Area].';
    areaTri(l) = max(Area);
    uh = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunction(geom, u, gradu, uh);
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






