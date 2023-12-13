clear all
close all
clc

%% Eseguo Triangolazione sul Dominio
area = 0.02;
geom = Triangolator(area);
close all

%% Problema differenziale
u = @(t,x,y) t + sin(x+y);
dtu= @(t,x,y) 1;
gradu = @(t,x,y) [cos(x+y), cos(x+y)]';
d2ux = @(t,x,y) [-sin(x + y), -sin(x + y)];
d2uy = @(t,x,y) [-sin(x + y), -sin(x + y)];
Hu = @(t,x,y) [d2ux(t,x,y); d2uy(t,x,y)]';
d2u = @(t,x,y) [1,0]*Hu(t,x,y)*[1,0]'+ [0,1]*Hu(t,x,y)*[0,1]';
rho = @(x,y) 1;
mu = @(x,y) 1.0;
beta = @(x,y) [3.0, 0.0];
sigma = @(x,y) 4.0;
f = @(t,x,y) dtu(t,x,y) - mu(x,y)*d2u(t,x,y)+beta(x,y)*gradu(t,x,y)+sigma(x,y)*u(t,x,y);
% f = @(x,y) 32*(x*(1-x) + y*(1-y));% + 16*(1-2*x)*(y*(1-y));
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(t,x,y) mu(x,y)*(n'*gradu(t,x,0));% @(x,y) -16*x*(1-x);
% gNe = @(x,y) -16*x*(1-x);
gDi = @(t,x,y) t + sin(x+y);
dtgDi = @(t,x,y) 1;
u0 = @(x,y) u(0,x,y);

%% Soluzione problema discretizzato
deltat = 0.01;
Nt = 5;
T = Nt*deltat;
uh = ParabolicP1(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);

%% Plot soluzione approssimata e esatta
% tTable = tTableforP2plot(geom.elements.triangles);
XY = geom.elements.coordinates;
x = XY(:,1);
y = XY(:,2);
Np = length(x);
soluzioneEsatta = zeros(Np,1);
for i=1:Np
    soluzioneEsatta(i) = u(T,x(i), y(i));
end
figure(1)
tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
trisurf(tTable, x, y, uh(:, Nt+1));
title("Grafico funzione approssimata")
figure(2)
trisurf(tTable, x, y, soluzioneEsatta);
title("Grafico soluzione esatta")
