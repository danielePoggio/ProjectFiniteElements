clear all
close all
clc

%% Problema differenziale
f = @(x,a) 1./(1 + exp(-a*(x-0.5)));
u = @(x,y) f(x+y,10);
% u = @(x,y) 16*x*(1-x)*y*(1-y);
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\calculateDerivate.m")
gradu = @(x,y) gradu(x,y)';
d2u = @(x,y) [1,0]*Hu(x,y)*[1,0]'+ [0,1]*Hu(x,y)*[0,1]';
mu = @(x,y) 0.01;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 500.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,0);
gDi = @(x,y) u(x,y);

%% Soluzione problema discretizzato
% Pk = 1;
% area = 0.002;
% geom = TriangolatorP2(area, Pk);
% close all
% % Siccome i coeff. sono costanti valutiamo Pe in maniera approssimata:
% Area = [geom.support.TInfo.Area].';
% area = max(Area);
% Se = sigma(0,0)*2*area/(6*mu(0,0));
% disp(Se)
% 
% % calcoliamo soluzione approssimata
% % uh = massLumping(geom, mu, sigma, f, gDi, gNe);
% uh = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
% XY = geom.elements.coordinates;
% ele = geom.elements.triangles;
% x = XY(:,1);
% y = XY(:,2);
% figure(1)
% tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
% % tTable = tTableforP2plot(ele);
% trisurf(tTable, x, y, uh);
% title("Grafico funzione approssimata")

%% Plot soluzione esatta
% sol = zeros(length(x),1);
% for i = 1:length(x)
%     sol(i) = u(x(i),y(i));
% end
% figure(2);
% trisurf(tTable, x, y, sol);
% xlabel('X');
% ylabel('Y');
% zlabel('u(x,y)');
% title("Soluzione esatta u(x,y)")

%% Valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
Pk = 1;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorL2vec = zeros(Ktest,1);
errorH1vec = zeros(Ktest,1);
Sevec = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end
    geom = TriangolatorP2(area,Pk);
    Area = [geom.support.TInfo.Area].';
    areaTri(l) = max(Area);
    h = sqrt(areaTri(l));
    mk = 1/24;
    Se = sigma(0,0)*2*area/(6*mu(0,0));
    Sevec(l) = Se;
    uh = massLumping(geom, mu, sigma, f, gDi, gNe);
%     uh = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
    % facciamo plot temporaneo della soluzione
    XY = geom.elements.coordinates;
    ele = geom.elements.triangles;
    x = XY(:,1);
    y = XY(:,2);
    figure(2+l)
    tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
    trisurf(tTable, x, y, uh);
    title("Grafico funzione approssimata")

    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
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