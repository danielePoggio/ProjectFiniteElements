clear all
close all 
clc

%% Condizioni di Dirichlet nulle
% definizione del problema e della soluzione
u = @(x,y) 16*x*(1-x)*y*(1-y);
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\calculateDerivate.m")
gradu = @(x,y) gradu(x,y)';
d2u = @(x,y) [1,0]*Hu(x,y)*[1,0]'+ [0,1]*Hu(x,y)*[0,1]';
mu = @(x,y) 1;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,0);
gDi = @(x,y) u(x,y);
% ordine di convergenza
Pk = 2;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorP2L2DiN = zeros(Ktest,1);
errorP2H1DiN = zeros(Ktest,1);
condAP2DiN = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end
    geom = TriangolatorP2Di(area);
    close all
    Area = [geom.support.TInfo.Area].';
    areaTri(l) = max(Area);
    [uh, condA] = FEMDiNeP2(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorP2L2DiN(l) = errorL2;
    errorP2H1DiN(l) = errorH1;
    condAP2DiN(l) = condA;
end
figure(1)
loglog(sqrt(areaTri), errorP2L2DiN)
title("Andamento errore norma L2")
saveas(1, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\L2P2DiNulle', "png");

pL2DiN = polyfit(log(sqrt(areaTri)), log(errorP2L2DiN), 1);

figure(2)
loglog(sqrt(areaTri), errorP2H1DiN)
title("Andamento errore norma H1")
saveas(2, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\H1P2DiNulle', "png");

pH1DiN = polyfit(log(sqrt(areaTri)), log(errorP2H1DiN), 1);

figure(3)
loglog(sqrt(areaTri), condAP2DiN)
title("Andamento condizionamento matrice di rigidezza")
saveas(3, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\CondAP2DiNulle', "png");

pcondADiN = polyfit(log(sqrt(areaTri)), log(condAP2DiN), 1);

save("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\P2DiNulle.mat","areaTri","errorP2L2DiN","errorP2H1DiN", "pL2DiN", "pH1DiN", "condAP2DiN", "pcondADiN" )

% Plot soluzione approssimata
% XY = geom.elements.coordinates;
% ele = geom.elements.triangles;
% x = XY(:,1);
% y = XY(:,2);
% figure(1)
% tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
% % tTable = tTableforP2plot(ele);
% trisurf(tTable, x, y, uh);
% title("Grafico funzione approssimata")
%% Condizioni di Dirichlet al bordo
% definizione problema differenziale
u = @(x,y) 16*x*(1-x)*y*(1-y)+x+y;
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\calculateDerivate.m")
gradu = @(x,y) gradu(x,y)';
d2u = @(x,y) [1,0]*Hu(x,y)*[1,0]'+ [0,1]*Hu(x,y)*[0,1]';
mu = @(x,y) 1.0;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*(n'*gradu(x,y));
gDi = @(x,y) u(x,y);
% ordine di convergenza
%% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
Pk = 2;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorP2L2DiNN = zeros(Ktest,1);
errorP2H1DiNN = zeros(Ktest,1);
condAP2DiNN = zeros(Ktest, 1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end
    clear geom
    geom = TriangolatorP2Di(area);
    close all
    Area = [geom.support.TInfo.Area].';
    maxArea = max(Area);
    areaTri(l) = maxArea;
    [uh, condA] = FEMDiNeP2(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorP2L2DiNN(l) = errorL2;
    errorP2H1DiNN(l) = errorH1; 
    condAP2DiNN(l) = condA;
end
figure(1)
loglog(sqrt(areaTri), errorP2L2DiNN)
title("Andamento errore norma L2")
saveas(1, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\L2P2DiNN', "png");

pL2DiNN = polyfit(log(sqrt(areaTri)), log(errorP2L2DiNN), 1);

figure(2)
loglog(sqrt(areaTri), errorP2H1DiNN)
title("Andamento errore norma H1")
saveas(2, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\H1P2DiNN', "png");

pH1DiNN = polyfit(log(sqrt(areaTri)), log(errorP2H1DiNN), 1);

figure(3)
loglog(sqrt(areaTri), condAP2DiNN)
title("Andamento condizionamento matrice di rigidezza")
saveas(3, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\CondAP2DiNN', "png");

pcondADiNN = polyfit(log(sqrt(areaTri)), log(condAP2DiN), 1);

save("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\P2DiNulle.mat","areaTri","errorP2L2DiN","errorP2H1DiN", "pL2DiN", "pH1DiN", "condAP2DiNN", "pcondADiNN" )

% Plot soluzione approssimata
% XY = geom.elements.coordinates;
% ele = geom.elements.triangles;
% x = XY(:,1);
% y = XY(:,2);
% figure(3)
% tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
% % tTable = tTableforP2plot(ele);
% trisurf(tTable, x, y, uh);
% 
% title("Grafico funzione approssimata")


%% Condizioni al bordo di Neumann
% definizione problema differenziale
u = @(x,y) 16*x*(1-x)*y*(1-y);
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\calculateDerivate.m")
gradu = @(x,y) gradu(x,y)';
d2u = @(x,y) [1,0]*Hu(x,y)*[1,0]'+ [0,1]*Hu(x,y)*[0,1]';
mu = @(x,y) 1;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,0);
gDi = @(x,y) u(x,y);
% ordine di convergenza
% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
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
loglog(sqrt(areaTri), errorP2L2DiNe)
title("Andamento errore norma L2")
saveas(1, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\L2P2DiNe', "png");

pL2DiNe = polyfit(log(sqrt(areaTri)), log(errorP2L2DiNe), 1);

figure(2)
loglog(sqrt(areaTri), errorP2H1DiNe)
title("Andamento errore norma H1")
saveas(2, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\H1P2DiNe', "png");

pH1DiNe = polyfit(log(sqrt(areaTri)), log(errorP2H1DiNe), 1);

figure(3)
loglog(sqrt(areaTri), condAP2DiNe)
title("Andamento condizionamento matrice di rigidezza")
saveas(3, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\condAP2DiNe', "png");

pcondADiNe = polyfit(log(sqrt(areaTri)), log(condAP2DiNe), 1);

save("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\P2DiNN.mat","areaTri","errorP2L2DiNe","errorP2H1DiNe", "condAP2DiNe", "pL2DiNe", "pH1DiNe", "pcondADiNe" )


% Plot soluzione approssimata
XY = geom.elements.coordinates;
ele = geom.elements.triangles;
x = XY(:,1);
y = XY(:,2);
figure(3)
tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
% tTable = tTableforP2plot(ele);
trisurf(tTable, x, y, uh);

title("Grafico funzione approssimata")


