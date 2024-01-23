clear all
close all
clc

%% Eseguo Triangolazione sul Dominio
area = 0.002;
geom = TriangolatorP1Di(area);
geom = P3(geom);

%% Problema differenziale
u = @(x,y) 16*x*(1-x)*y*(1-y);
[gradu, d2u] = calculateDerivate(u);
mu = @(x,y) 1;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
% f = @(x,y) 32*(x*(1-x) + y*(1-y));
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,y);
% gNe = @(x,y) -16*x*(1-x)-1;
gDi = @(x,y) u(x,y);

%% Assemblaggio
[uh, idxValue] = FEMDiNeP3(geom, mu, beta, sigma, f, gDi, gNe);

%% estrazione valori ai vertici:
uhpoint = uh(idxValue);
XY = geom.elements.coordinates;
ele = geom.elements.triangles;
x = XY(:,1);
y = XY(:,2);
figure(1)
tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
trisurf(tTable, x, y, uhpoint);
title("Grafico funzione approssimata")

