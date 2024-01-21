clear all
close all
clc

%% Definizione problema differenziale
u = @(t,x,y) x+y+sin(t);
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\calculateDerivateTemporal.m")
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

%% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
Pk = 1;
area = 0.01;
geom = TriangolatorP1Di(area);
close all
Area = [geom.support.TInfo.Area].';
h = sqrt(max(Area));

%%  Eseguo Test Crank Nicolson
XY = geom.elements.coordinates;
x = XY(:,1);
y = XY(:,2);
Np = length(x);
T = 2.0;
Ktest = 4;
deltaTest = zeros(Ktest,1);
deltaTest(1) = 0.1;
errorL2CN = zeros(Ktest,1);
errorH1CN = zeros(Ktest,1);
errorLInfCN = zeros(Ktest,1);
condBCN = zeros(Ktest,1);
numberStep = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        deltat = deltaTest(1);
    else
        deltat = deltaTest(l-1)/4;
        deltaTest(l) = deltat;
    end
    Nt = T/deltat;
    uT = @(x,y) u(T,x,y);
    graduT = @(x,y) gradu(T,x,y);
    soluzioneEsatta = zeros(Np,1);
    for i=1:Np
        soluzioneEsatta(i) = u(T,x(i), y(i));
    end
    numberStep(l) = Nt;
    [uh, condB] = crankNicolson(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);
    condBCN(l) = condB;
    uhT = uh(:,Nt+1);
    [errorL2, errorH1] = errorFunctionOld(geom, uT, graduT, uhT, Pk);
    errorL2CN(l) = errorL2;
    errorH1CN(l) = errorH1;
    errorLInfCN(l) = norm(soluzioneEsatta - uhT, 'inf');
end

% figure(1)
% loglog(deltaTest, errorH1CN)
% title("Andamento errore norma H1")

pH1CN = polyfit(log(deltaTest), log(errorH1CN), 1);

% figure(2)
% loglog(deltaTest, errorL2CN)
% title("Andamento errore norma L2")

pL2CN = polyfit(log(deltaTest), log(errorL2CN), 1);

condBCN = polyfit(log(deltaTest), log(condBCN), 1);

% figure(3)
% plot(deltaTest, errorLInfCN)
% title("Andamento errore norma LInf")

%%  Eseguo Test Eulero Implicito
XY = geom.elements.coordinates;
x = XY(:,1);
y = XY(:,2);
Np = length(x);
T = 2.0;
Ktest = 4;
deltaTest = zeros(Ktest,1);
deltaTest(1) = 0.1;
errorL2IE = zeros(Ktest,1);
errorH1IE = zeros(Ktest,1);
condBIE = zeros(Ktest,1);
errorLInfIE = zeros(Ktest,1);
numberStep = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        deltat = deltaTest(1);
    else
        deltat = deltaTest(l-1)/4;
        deltaTest(l) = deltat;
    end
    Nt = T/deltat;
    uT = @(x,y) u(T,x,y);
    graduT = @(x,y) gradu(T,x,y);
    soluzioneEsatta = zeros(Np,1);
    for i=1:Np
        soluzioneEsatta(i) = u(T,x(i), y(i));
    end
    numberStep(l) = Nt;
    [uh, condB] = implicitEuler(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);
    condBIE(l) = condB;
    uhT = uh(:,Nt+1);
    [errorL2, errorH1] = errorFunctionOld(geom, uT, graduT, uhT, Pk);
    errorL2IE(l) = errorL2;
    errorH1IE(l) = errorH1;
    errorLInfIE(l) = norm(soluzioneEsatta - uhT, 'inf');
end

figure(1)
loglog(deltaTest, errorH1IE)
hold on
loglog(deltaTest, errorH1CN)
title("Andamento errore norma H1")
legend("Errore Eulero Implicito", "Errore Crank-Nicolson")
xlabel("delta t")
ylabel("Errore in norma H1")

pH1IE = polyfit(log(deltaTest), log(errorH1IE), 1);

figure(2)
loglog(deltaTest, errorL2IE)
hold on
loglog(deltaTest, errorL2CN)
title("Andamento errore norma L2 IE")
legend("Errore Eulero Implicito", "Errore Crank-Nicolson")
xlabel("x")
ylabel("Errore in norma L2")

pL2IE = polyfit(log(deltaTest), log(errorL2IE), 1);

condBIE = polyfit(log(deltaTest), log(condBIE), 1);


% figure(Ktest+4)
% plot(deltaTest, errorLInfIE)
% title("Andamento errore norma LInf")