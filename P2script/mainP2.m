clear all
close all
clc
% ciao dfhfgd
%% Eseguo Triangolazione sul Dominio
area = 0.002;
geom = Triangolator(area);
run("P2")
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
uh = FEMDiNeP2(geom, mu, beta, sigma, f, gDi, gNe);

%% Plot soluzione approssimata
XY = geom.elements.coordinates;
x = XY(:,1);
y = XY(:,2);
figure(1)
tri = delaunay(x, y); % Genera la matrice di connettività dei triangoli
trisurf(tri, x, y, uh);
title("Grafico funzione approssimata")

