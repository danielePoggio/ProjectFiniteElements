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
Pk = 1;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorL2DiN = zeros(Ktest,1);
errorH1DiN = zeros(Ktest,1);
condADiN = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end
    geom = TriangolatorP1Di(area);
    close all
    Area = [geom.support.TInfo.Area].';
    areaTri(l) = max(Area);
    [uh, condA] = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorL2DiN(l) = errorL2;
    errorH1DiN(l) = errorH1;
    condADiN(l) = condA;
end
figure(1)
loglog(sqrt(areaTri), errorL2DiN)
title("Andamento errore norma L2")
saveas(1, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\L2DiNulle', "png");

pL2DiN = polyfit(log(sqrt(areaTri)), log(errorL2DiN), 1);

figure(2)
loglog(sqrt(areaTri), errorH1DiN)
title("Andamento errore norma H1")
saveas(2, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\H1DiNulle', "png");

pH1DiN = polyfit(log(sqrt(areaTri)), log(errorH1DiN), 1);

figure(3)
loglog(sqrt(areaTri), condADiN)
title("Andamento condizionamento matrice di rigidezza")
saveas(3, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\CondADiNulle', "png");

pcondADiN = polyfit(log(sqrt(areaTri)), log(condADiN), 1);

save("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\P1DiNulle.mat","areaTri","errorL2DiN","errorH1DiN", "pL2DiN", "pH1DiN", "condADiN", "pcondADiN" )

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
Pk = 1;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorL2DiNN = zeros(Ktest,1);
errorH1DiNN = zeros(Ktest,1);
condADiNN = zeros(Ktest, 1);
for l=1:Ktest
    if l == 1
        area = areaTri(1);
    else
        area = areaTri(l-1)/4;
    end
    clear geom
    geom = TriangolatorP1Di(area);
    close all
    Area = [geom.support.TInfo.Area].';
    maxArea = max(Area);
    areaTri(l) = maxArea;
    [uh, condA] = FEMDiNe(geom, mu, beta, sigma, f, gDi, gNe);
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorL2DiNN(l) = errorL2;
    errorH1DiNN(l) = errorH1; 
    condADiNN(l) = condA;
end
figure(1)
loglog(sqrt(areaTri), errorL2DiNN)
title("Andamento errore norma L2")
saveas(1, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\L2DiNN', "png");

pL2DiNN = polyfit(log(sqrt(areaTri)), log(errorL2DiNN), 1);

figure(2)
loglog(sqrt(areaTri), errorH1DiNN)
title("Andamento errore norma H1")
saveas(2, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\H1DiNN', "png");

pH1DiNN = polyfit(log(sqrt(areaTri)), log(errorH1DiNN), 1);

figure(3)
loglog(sqrt(areaTri), condADiNN)
title("Andamento condizionamento matrice di rigidezza")
saveas(3, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\CondADiNN', "png");

pcondADiNN = polyfit(log(sqrt(areaTri)), log(condADiN), 1);

save("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\P1DiNulle.mat","areaTri","errorL2DiN","errorH1DiN", "pL2DiN", "pH1DiN", "condADiNN", "pcondADiNN" )

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
u = @(x,y) 16*x*(1-x)*y*(1-y)+x+y;
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\calculateDerivate.m")
gradu = @(x,y) gradu(x,y)';
d2u = @(x,y) [1,0]*Hu(x,y)*[1,0]'+ [0,1]*Hu(x,y)*[0,1]';
mu = @(x,y) 1;
beta = @(x,y) [x, 0.0];
sigma = @(x,y) -y;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*(n'*gradu(x,y));
gDi = @(x,y) u(x,y);
% ordine di convergenza
% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
Pk = 1;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorL2DiNe = zeros(Ktest,1);
errorH1DiNe = zeros(Ktest,1);
condAvec = zeros(Ktest,1);
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
    condAvec(l) = condA;
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorL2DiNe(l) = errorL2;
    errorH1DiNe(l) = errorH1; 
end
figure(1)
loglog(sqrt(areaTri), errorL2DiNe)
title("Andamento errore norma L2")
saveas(1, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\L2DiNe', "png");

pL2DiNe = polyfit(log(sqrt(areaTri)), log(errorL2DiNe), 1);

figure(2)
loglog(sqrt(areaTri), errorH1DiNe)
title("Andamento errore norma H1")
saveas(2, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\H1DiNe', "png");

pH1DiNe = polyfit(log(sqrt(areaTri)), log(errorH1DiNe), 1);

figure(3)
loglog(sqrt(areaTri), condAvec)
title("Andamento condizionamento matrice di rigidezza")
saveas(3, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\condADiNe', "png");

pcondADiNe = polyfit(log(sqrt(areaTri)), log(condAvec), 1);

save("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\P1DiNN.mat","areaTri","errorL2DiNe","errorH1DiNe", "condAvec", "pL2DiNe", "pH1DiNe", "pcondADiNe" )


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

