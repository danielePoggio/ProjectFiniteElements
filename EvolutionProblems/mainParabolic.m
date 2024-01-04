clear all
close all
clc

%% Eseguo Triangolazione sul Dominio
% area = 0.02;
% geom = Triangolator(area);
% close all

%% Problema differenziale
u = @(t,x,y) x+y+t^5;
% gradu = @(t,x,y) [cos(x+y), cos(x+y)]';
% u = @(t,x,y) x+y + sin(t^2);
run("calculateDerivate.m")
d2u = @(t,x,y) [1,0]*Hu(t,x,y)*[1,0]'+ [0,1]*Hu(t,x,y)*[0,1]';
rho = @(x,y) 1;
mu = @(x,y) 1.0;
beta = @(x,y) [3.0, 0.0];
sigma = @(x,y) 4.0;
f = @(t,x,y) dut(t,x,y) - mu(x,y)*d2u(t,x,y)+beta(x,y)*gradu(t,x,y)+sigma(x,y)*u(t,x,y);
% f = @(x,y) 32*(x*(1-x) + y*(1-y));% + 16*(1-2*x)*(y*(1-y));
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(t,x,y) mu(x,y)*(n'*gradu(t,x,0));% @(x,y) -16*x*(1-x);
% gNe = @(x,y) -16*x*(1-x);
gDi = @(t,x,y) t^2 + sin(x+y);
dtgDi = @(t,x,y) 2*t;
u0 = @(x,y) u(0,x,y);

%% Soluzione problema discretizzato
% deltat = 0.01;
% Nt = 10;
% T = Nt*deltat;
% uh = ParabolicP1(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);

%% Plot soluzione approssimata e esatta
% tTable = tTableforP2plot(geom.elements.triangles);
XY = geom.elements.coordinates;
x = XY(:,1);
y = XY(:,2);
Np = length(x);
soluzioneEsatta = zeros(Np,1);
for i=1:Np
    soluzioneEsatta(i) = u(T,x(i), y(i));
end
figure(1)
tTable = delaunay(x, y); % Genera la matrice di connettività dei triangoli
trisurf(tTable, x, y, uh(:, Nt+1));
title("Grafico funzione approssimata")
figure(2)
trisurf(tTable, x, y, soluzioneEsatta);
title("Grafico soluzione esatta")

%% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
% TEST SUL PASSO TEMPORALE
Pk = 1;
T = 0.5;
% calcolo problema stazionario all'istante finale T
uT = @(x,y) u(T,x,y);
graduT = @(x,y) gradu(T,x,y);
HuT = @(x,y) Hu(T,x,y);
d2uT = @(x,y) d2u(T,x,y);
gNeT = @(x,y) gNe(T,x,y);% @(x,y) -16*x*(1-x);
% gNe = @(x,y) -16*x*(1-x);
gDiT = @(x,y) gDi(T,x,y);
dtgDiT = @(x,y) dtgDi(T,x,y);
fT = @(x,y) f(T,x,y);

% calcoliamo ora errore rispetto al problema parabolico
clear geom
area = 0.005;
geom = Triangolator(area);
close all
Area = [geom.support.TInfo.Area].';
h = sqrt(max(Area));
Ktest = 5;
deltaTest = zeros(Ktest,1);
deltaTest(1) = 0.01;
errorC0vec = zeros(Ktest,1);
errorL2timevec = zeros(Ktest,1);
errorL2vec = zeros(Ktest,1);
for l=1:Ktest
    if l == 1
        deltat = deltaTest(1);
    else
        deltat = deltaTest(l-1)/2;
        deltaTest(l) = deltat;
    end
    Nt = T/deltat;
    uh = ParabolicP1(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0);
    errorC0 = 0;
    errorL2timeSquared = 0;
    for n=1:Nt
        T_ = n*deltat;
        ut = @(x,y) u(T_,x,y);
        gradut = @(x,y) gradu(T_,x,y);
        [errorL2, errorH1] = errorFunction(geom, ut, gradut, uh(:, n), Pk);
        if errorC0 < errorL2
            errorC0 = errorL2;
        end
        errorL2timeSquared = errorL2timeSquared + (errorH1^2)*deltat;
        if n == Nt
            errorL2vec(l) = errorL2;
        end

    end
    errorC0vec(l) = errorC0;
    errorL2timevec(l) = sqrt(errorL2timeSquared);
end

figure(1)
plot(deltaTest, errorC0vec+errorL2timevec)
title("Andamento errore completo")

figure(2)
plot(deltaTest, errorL2vec)
title("Andamento errore norma L2")

pL2 = polyfit(log(deltaTest), log(errorC0vec), 1);

% figure(2)
% loglog(deltaTest, errorH1vec-errorT(2))
% title("Andamento errore norma H1")
% 
% pH1 = polyfit(log(deltaTest), log(errorL2timevec), 1);

%% Plot soluzione ai vari istanti di tempo
% for j=1:Nt+1
%     figure(j)
%     XY = geom.elements.coordinates;
%     x = XY(:,1);
%     y = XY(:,2);
%     tTable = delaunay(x, y);
%     trisurf(tTable, x, y, uh(:, j));
% end