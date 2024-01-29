clear all
close all
clc

%% Confronto P1 e P2 su stessa soluzione:
u = @(x,y) 16*x*(1-x)*y*(1-y);
[gradu, d2u] = calculateDerivate(u);
mu = @(x,y) 1;
beta = @(x,y) [3.0, 1.0];
sigma = @(x,y) 2;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,0);
gDi = @(x,y) u(x,y);

%% P1:
Pk = 1;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorP1L2DiNe = zeros(Ktest,1);
errorP1H1DiNe = zeros(Ktest,1);
condAP1DiNe = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end
    clear geom
    geom = TriangolatorP1Ne(area);
    close all
    Area = [geom.support.TInfo.Area].';
    maxArea = max(Area);
    areaTri(l) = maxArea;
    [uh, condA] = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
    condAP1DiNe(l) = condA;
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorP1L2DiNe(l) = errorL2;
    errorP1H1DiNe(l) = errorH1; 
end

%% P2
Pk = 2;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorP2L2DiNe = zeros(Ktest,1);
errorP2H1DiNe = zeros(Ktest,1);
condAP2DiNe = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end
    clear geom
    geom = TriangolatorP2Ne(area);
    close all
    Area = [geom.support.TInfo.Area].';
    maxArea = max(Area);
    areaTri(l) = maxArea;
    [uh, condA] = FEMDiNeP2(geom, mu, beta, sigma, f, gDi, gNe);
    condAP2DiNe(l) = condA;
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorP2L2DiNe(l) = errorL2;
    errorP2H1DiNe(l) = errorH1; 
end

figure(1)
loglog(sqrt(areaTri), errorP1L2DiNe)
hold on
loglog(sqrt(areaTri), errorP2L2DiNe)
legend('Errore norma L2 P1', 'Errore norma L2 P2');
xlabel('h')
ylabel('Errore norma L2')
hold off

figure(2)
loglog(sqrt(areaTri), errorP1H1DiNe)
hold on
loglog(sqrt(areaTri), errorP2H1DiNe)
legend('Errore norma H1 P1', 'Errore norma H1 P2');
xlabel('h')
ylabel('Errore norma H1')
hold off

figure(3)
loglog(sqrt(areaTri), condAP1DiNe)
hold on
loglog(sqrt(areaTri), condAP2DiNe)
legend('Condizionamento A per P1', 'Condizionamento A per P2');
xlabel('h')
ylabel('Numero di condizionamento')
hold off