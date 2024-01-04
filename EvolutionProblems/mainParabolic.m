clear all
close all
clc

%% Eseguo Triangolazione sul Dominio
% area = 0.02;
% geom = Triangolator(area);
% close all

%% Problema differenziale
u = @(t,x,y) x+y+sin(t);
run("calculateDerivateTemporal.m")
gradu = @(t,x,y) gradu(t,x,y)';
d2u = @(t,x,y) [1,0]*Hu(t,x,y)*[1,0]'+ [0,1]*Hu(t,x,y)*[0,1]';
rho = @(x,y) 1.0;
mu = @(x,y) 1.0;
beta = @(x,y) [3.0, 0.0];
sigma = @(x,y) 4.0;
f = @(t,x,y) dut(t,x,y)-mu(x,y)*d2u(t,x,y)+beta(x,y)*gradu(t,x,y)+sigma(x,y)*u(t,x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(t,x,y) mu(x,y)*(n'*gradu(t,x,0));% @(x,y) -16*x*(1-x);
% gNe = @(x,y) -16*x*(1-x);
gDi = @(t,x,y) u(t,x,y);
dtgDi = @(t,x,y) dut(t,x,y);
u0 = @(x,y) u(0,x,y);

%% Soluzione problema discretizzato
% deltat = 0.01;
% Nt = 10;
% T = Nt*deltat;
% uh = ParabolicP1new(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);

%% Plot soluzione approssimata e esatta
% tTable = tTableforP2plot(geom.elements.triangles);
% XY = geom.elements.coordinates;
% x = XY(:,1);
% y = XY(:,2);
% Np = length(x);
% soluzioneEsatta = zeros(Np,1);
% for i=1:Np
%     soluzioneEsatta(i) = u(T,x(i), y(i));
% end
% figure(1)
% tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
% trisurf(tTable, x, y, uh(:, Nt+1));
% title("Grafico funzione approssimata")
% figure(2)
% trisurf(tTable, x, y, soluzioneEsatta);
% title("Grafico soluzione esatta")

%% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
% TEST SUL PASSO TEMPORALE
Pk = 1;
T = 2.0;
uT = @(x,y) u(T,x,y);
graduT = @(x,y) gradu(T,x,y);

% calcoliamo ora errore rispetto al problema parabolico
clear geom
area = 0.02;
geom = Triangolator(area);
close all
Area = [geom.support.TInfo.Area].';
h = sqrt(max(Area));
Ktest = 3;
deltaTest = zeros(Ktest,1);
deltaTest(1) = 0.01;
errorL2vec = zeros(Ktest,1);
errorH1vec = zeros(Ktest,1);
numberStep = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        deltat = deltaTest(1);
    else
        deltat = deltaTest(l-1)/4;
        deltaTest(l) = deltat;
    end
    Nt = T/deltat;
    numberStep(l) = Nt;
    uh = ParabolicP1new(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);
    [errorL2, errorH1] = errorFunctionNew(geom, uT, graduT, uh(:,Nt+1), Pk);
    errorL2vec(l) = errorL2;
    errorH1vec(l) = errorH1;
end

figure(1)
plot(deltaTest, errorH1vec)
title("Andamento errore norma H1")

figure(2)
plot(deltaTest, errorL2vec)
title("Andamento errore norma L2")



%% Plot soluzione ai vari istanti di tempo
for j=1:Nt+1
    figure(j)
    XY = geom.elements.coordinates;
    x = XY(:,1);
    y = XY(:,2);
    tTable = delaunay(x, y);
    trisurf(tTable, x, y, uh(:, j));
end

XY = geom.elements.coordinates;
x = XY(:,1);
y = XY(:,2);
Np = length(x);
soluzioneEsatta = zeros(Np,1);
for i=1:Np
    soluzioneEsatta(i) = u(T,x(i), y(i));
end
figure(Nt+2)
trisurf(tTable, x, y, soluzioneEsatta);
title("Grafico soluzione esatta")