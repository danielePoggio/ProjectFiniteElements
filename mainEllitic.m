clear all
close all
clc

%% Eseguo Triangolazione sul Dominio
area = 0.02;
geom = Triangolator(area);
close all

%% Problema differenziale
u = @(x,y) 16*x*(1-x)*y*(1-y);
gradu = @(x,y) 16*[y*(1-y)*(1 -2*x), x*(1-x)*(1 - 2*y)]';
d2ux = @(x,y) [32*conj(y)*(conj(y) - 1), 16*(conj(x) - 1)*(2*conj(y) - 1) + 16*conj(x)*(2*conj(y) - 1)];
d2uy = @(x,y) [16*(conj(y) - 1)*(2*conj(x) - 1) + 16*conj(y)*(2*conj(x) - 1),32*conj(x)*(conj(x) - 1)];
Hu = @(x,y) [d2ux(x,y); d2uy(x,y)]';
d2u = @(x,y) [1,0]*Hu(x,y)*[1,0]'+ [0,1]*Hu(x,y)*[0,1]';
mu = @(x,y) 0.001;
beta = @(x,y) [3.0, 0.0];
sigma = @(x,y) 0.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
% f = @(x,y) 32*(x*(1-x) + y*(1-y));% + 16*(1-2*x)*(y*(1-y));
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) n'*gradu(x,0);% @(x,y) -16*x*(1-x);
% gNe = @(x,y) -16*x*(1-x);
gDi = @(x,y) 0;


%% Soluzione problema discretizzato
% uh = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
Pk = 1;
uh = SUPG(geom, Pk, mu, beta, f, gDi, gNe);
%% Plot soluzione approssimata
XY = geom.elements.coordinates;
ele = geom.elements.triangles;
x = XY(:,1);
y = XY(:,2);
figure(1)
tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
% tTable = tTableforP2plot(ele);
trisurf(tTable, x, y, uh);
title("Grafico funzione approssimata")
% 
%% calcoliamo stima dell'errore
% [errorL2, errorH1] = errorFunction(geom, u, gradu, uh, 1);
% 
% %% Andiamo a vedere come estrarre il valore minimo e massimo dell'area nella triangolazione
% Area = [geom.support.TInfo.Area].';
% areaMax = max(Area);

%% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
% Pk = 1;
% Ktest = 4;
% areaTri = zeros(Ktest,1);
% areaTri(1) = 0.01;
% errorL2vec = zeros(Ktest,1);
% errorH1vec = zeros(Ktest,1);
% for l=1:Ktest
%     if l == 1
%         area = areaTri(1);
%     else
%         area = areaTri(l-1)/4;
%     end
%     geom = Triangolator(area);
%     close all
%     Area = [geom.support.TInfo.Area].';
%     areaTri(l) = max(Area);
%     uh = FEMDiNeQuadratura(geom, mu, beta, sigma, f, gDi, gNe);
%     [errorL2, errorH1] = errorFunction(geom, u, gradu, uh, Pk);
%     errorL2vec(l) = errorL2;
%     errorH1vec(l) = errorH1;
% end
% figure(1)
% loglog(sqrt(areaTri), errorL2vec)
% title("Andamento errore norma L2")
% 
% pL2 = polyfit(log(sqrt(areaTri)), log(errorL2vec), 1);
% 
% figure(2)
% loglog(sqrt(areaTri), errorH1vec)
% title("Andamento errore norma H1")
% 
% pH1 = polyfit(log(sqrt(areaTri)), log(errorH1vec), 1);
