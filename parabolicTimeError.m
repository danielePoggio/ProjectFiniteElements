clear all
close all
clc

%% Problema differenziale
u = @(t,x,y) x+y+sin(t);
run("calculateDerivateTemporal.m")
d2u = @(t,x,y) [1,0]*Hu(t,x,y)*[1,0]'+ [0,1]*Hu(t,x,y)*[0,1]';
rho = @(x,y) 1.0;
mu = @(x,y) 1.0;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(t,x,y) rho(x,y)*dut(t,x,y)-mu(x,y)*d2u(t,x,y)+beta(x,y)*gradu(t,x,y)'+sigma(x,y)*u(t,x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(t,x,y) mu(x,y)*(n'*gradu(t,x,0)');
gDi = @(t,x,y) u(t,x,y);
dtgDi = @(t,x,y) dut(t,x,y);
u0 = @(x,y) u(0,x,y);

%% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
% TEST SUL PASSO TEMPORALE
Pk = 1;
T = 100.0;
uT = @(x,y) u(T,x,y);
graduT = @(x,y) gradu(T,x,y)';

%calcoliamo ora errore rispetto al problema parabolico
area = 0.002;
geom = Triangolator(area);
close all
Area = [geom.support.TInfo.Area].';
h = sqrt(max(Area));
% plot soluzione finale
XY = geom.elements.coordinates;
x = XY(:,1);
y = XY(:,2);
Np = length(x);
soluzioneEsatta = zeros(Np,1);
for i=1:Np
    soluzioneEsatta(i) = u(T,x(i), y(i));
end
tTable = delaunay(x, y);
figure(1)
trisurf(tTable, x, y, soluzioneEsatta);
title("Grafico soluzione esatta")

% Eseguo Test
Ktest = 3;
deltaTest = zeros(Ktest,1);
deltaTest(1) = 10;
errorL2vec = zeros(Ktest,1);
errorH1vec = zeros(Ktest,1);
errorLInfvec = zeros(Ktest,1);
numberStep = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        deltat = deltaTest(1);
    else
        deltat = deltaTest(l-1)/2;
        deltaTest(l) = deltat;
    end
    Nt = T/deltat;
    numberStep(l) = Nt;
    uh = parabolicCN(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);
    uhT = uh(:,Nt+1);
    [errorL2, errorH1] = errorFunctionOld(geom, uT, graduT, uhT, Pk);
    errorL2vec(l) = errorL2;
    errorH1vec(l) = errorH1;
    errorLInfvec(l) = norm(soluzioneEsatta - uhT, 'inf');
    % plotto soluzione approssimata
    figure(l+1)
    trisurf(tTable, x, y, uh(:, Nt+1));
    title("Grafico soluzione approssimata")

end

figure(Ktest+2)
plot(deltaTest, errorH1vec)
title("Andamento errore norma H1")

pH1 = polyfit(log(deltaTest), log(errorH1vec), 1);

figure(Ktest+3)
plot(deltaTest, errorL2vec)
title("Andamento errore norma L2")

pL2 = polyfit(log(deltaTest), log(errorL2vec), 1);

figure(Ktest+4)
plot(deltaTest, errorLInfvec)
title("Andamento errore norma LInf")

