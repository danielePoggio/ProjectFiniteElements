%% Formule di quadratura per coeff non costanti
% definizione problema differenziale
u = @(x,y) 16*x*(1-x)*y*(1-y)+x+y;
[gradu, d2u] = calculateDerivate(u);
mu = @(x,y) 1;
beta = @(x,y) [x, 0.0];
sigma = @(x,y) -y;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*(n'*gradu(x,y));
gDi = @(x,y) u(x,y);

%% ordine di convergenza
% valutiamo come cambiano gli errori in norma L2 ed H1 al variare dell'area massima della triangolazione
Pk = 1;
Ktest = 3;
areaTri = zeros(Ktest,1);
areaTri(1) = 0.01;
errorL2DiNeQ = zeros(Ktest,1);
errorH1DiNeQ = zeros(Ktest,1);
condADiNeQ = zeros(Ktest,1);
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
    [uh, condA] = FEMDiNeQuadratura(geom, mu, beta, sigma, f, gDi, gNe);
    condADiNeQ(l) = condA;
    [errorL2, errorH1] = errorFunctionOld(geom, u, gradu, uh, Pk);
    errorL2DiNeQ(l) = errorL2;
    errorH1DiNeQ(l) = errorH1; 
end
figure(1)
loglog(sqrt(areaTri), errorL2DiNeQ)
title("Andamento errore norma L2")
saveas(1, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\L2DiNeQ', "png");

pL2DiNeQ = polyfit(log(sqrt(areaTri)), log(errorL2DiNe), 1);

figure(2)
loglog(sqrt(areaTri), errorH1DiNeQ)
title("Andamento errore norma H1")
saveas(2, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\H1DiNeQ', "png");

pH1DiNeQ = polyfit(log(sqrt(areaTri)), log(errorH1DiNe), 1);

figure(3)
loglog(sqrt(areaTri), condADiNeQ)
title("Andamento condizionamento matrice di rigidezza")
saveas(3, 'C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\condADiNeQ', "png");

pcondADiNeQ = polyfit(log(sqrt(areaTri)), log(condAvec), 1);

figure(4)
plot(sqrt(areaTri), errorL2DiNeQ)
hold on
plot(sqrt(areaTri), errorL2DiNe)
% Aggiunta della legenda
legend('Errore formule di quadratura', 'Errore approssimazione formule P1'); % 'best' posiziona la legenda nella posizione migliore

% Aggiunta di etichette agli assi e titolo al grafico
xlabel('h');
ylabel('Errore in norma L2');
hold off


errorL2DiNeQ - errorL2DiNe

save("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\Experiment\P1DiNeQ.mat","areaTri","errorL2DiNeQ","errorH1DiNeQ", "condADiNeQ", "pL2DiNeQ", "pH1DiNeQ", "pcondADiNeQ" )