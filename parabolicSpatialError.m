clear all
close all
clc

%% Problema differenziale
u = @(t,x,y) sin(x+y)+ t;
run("calculateDerivateTemporal.m")
d2u = @(t,x,y) [1,0]*Hu(t,x,y)*[1,0]'+ [0,1]*Hu(t,x,y)*[0,1]';
rho = @(x,y) 1.0;
mu = @(x,y) 1.0;
beta = @(x,y) [3.0, 0.0];
sigma = @(x,y) 4.0;
f = @(t,x,y) rho(x,y)*dut(t,x,y)-mu(x,y)*d2u(t,x,y)+beta(x,y)*gradu(t,x,y)'+sigma(x,y)*u(t,x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(t,x,y) mu(x,y)*(n'*gradu(t,x,0)');
gDi = @(t,x,y) u(t,x,y);
dtgDi = @(t,x,y) dut(t,x,y);
u0 = @(x,y) u(0,x,y);
%% Errore Spaziale
Pk = 1;
T = 2.0;
deltat = 0.1;
Nt = T/deltat;
uT = @(x,y) u(T,x,y);
graduT = @(x,y) gradu(T,x,y)';
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
    % creazione mesh
    geom = Triangolator(area);
    Area = [geom.support.TInfo.Area].';
    areaTri(l) = max(Area);

    %risoluzione del problema parabolico
    uh = parabolicCN(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);
    [errorL2, errorH1] = errorFunction(geom, uT, graduT, uh(:,Nt+1), Pk);

    % plot soluzione approssimata e soluzione finale
    XY = geom.elements.coordinates;
    x = XY(:,1);
    y = XY(:,2);
    Np = length(x);
    soluzioneEsatta = zeros(Np,1);
    for i=1:Np
        soluzioneEsatta(i) = u(T,x(i), y(i));
    end
    tTable = delaunay(x, y);
    figure(2*l-1)
    trisurf(tTable, x, y, uh(:, Nt+1));
    title("Grafico soluzione approssimata")

    figure(2*l)
    trisurf(tTable, x, y, soluzioneEsatta);
    title("Grafico soluzione esatta")
    errorL2vec(l) = errorL2;
    errorH1vec(l) = errorH1;
end

figure(2*Ktest+1)
loglog(sqrt(areaTri), errorL2vec)
title("Andamento errore norma L2")

pL2 = polyfit(log(sqrt(areaTri)), log(errorL2vec), 1);

figure(2*Ktest+2)
loglog(sqrt(areaTri), errorH1vec)
title("Andamento errore norma H1")

pH1 = polyfit(log(sqrt(areaTri)), log(errorH1vec), 1);

