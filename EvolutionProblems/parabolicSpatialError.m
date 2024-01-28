clear all
close all
clc

%% Problema differenziale
u = @(t,x,y) sin(x+y) + t.^2;
run("calculateDerivateTemporal.m")
rho = @(x,y) 1.0;
mu = @(x,y) 1.0;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(t,x,y) rho(x,y)*dut(t,x,y)-mu(x,y)*d2u(t,x,y)+beta(x,y)*gradu(t,x,y)+sigma(x,y)*u(t,x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(t,x,y) mu(x,y)*(n'*gradu(t,x,y));
gDi = @(t,x,y) u(t,x,y);
dtgDi = @(t,x,y) dut(t,x,y);
u0 = @(x,y) u(0,x,y);

%% Errore Spaziale
Pk = 2; % ordine dei polinomi
T = 5.0;
deltat = 0.5;
Nt = T/deltat;
uT = @(x,y) u(T,x,y);
graduT = @(x,y) gradu(T,x,y);
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorL2space = zeros(Ktest,1);
errorH1space = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end

    % creazione mesh
    geom = TriangolatorP2Ne(area);
    Area = [geom.support.TInfo.Area].';
    areaTri(l) = max(Area);

    %risoluzione del problema parabolico
    uh = crankNicolsonP2(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);
    uhT = uh(:,Nt+1);
    [errorL2, errorH1] = errorFunctionOld(geom, uT, graduT, uhT, Pk);
    errorL2space(l) = errorL2;
    errorH1space(l) = errorH1;
end

figure(1)
loglog(sqrt(areaTri), errorL2space)
hold on
loglog(sqrt(areaTri), errorH1space)
legend("L2", "H1")
xlabel("h")
ylabel("Errori")
title("Andamento errori")

pL2 = polyfit(log(sqrt(areaTri)), log(errorL2space), 1);

figure(2)
loglog(sqrt(areaTri), errorH1space)
title("Andamento errore norma H1")

pH1 = polyfit(log(sqrt(areaTri)), log(errorH1space), 1);

