function uh = FEMDiNeP2tentativo(geom, mu, beta, sigma, f, gDi, gNe)
%Assemblaggio FEM Dirichlet non omogeneo, calcolo coeff con nodi di
%quadratura invece che approssimare con il baricentro
pivot = geom.pivot.pivot;
ele = geom.elements.triangles;
XY = geom.elements.coordinates;
Np = length(XY);
Ndof = max(pivot);
NDi = -min(pivot);
Nele = length(ele);
A = zeros(Ndof,Ndof);
Ad = zeros(Ndof, NDi);
b = zeros(Ndof,1);

% ricaviamo nodi di quadratura
run("nodes_weights.m")
% estraiamo funzioni P2
Nv = 6;
N1 = @(x,y) x;
N2 = @(x,y) y;
N3 = @(x,y) 1 - x - y;
phi1 = @(x,y) 2*N1(x,y)*(N1(x,y) - 0.5);
phi2 = @(x,y) 2*N2(x,y)*(N2(x,y) - 0.5);
phi3 = @(x,y) 2*N3(x,y)*(N3(x,y) - 0.5);
phi4 = @(x,y) 4*N3(x,y)*N1(x,y);
phi5 = @(x,y) 4*N1(x,y)*N2(x,y);
phi6 = @(x,y) 4*N2(x,y)*N3(x,y);
phi = @(x,y) [phi1(x,y), phi2(x,y), phi3(x,y), phi4(x,y), phi5(x,y), phi6(x,y)]';

Jphi = @(x,y) [4*x - 1, 0, 4*x + 4*y - 3, 4 - 4*y - 8*x, 4*y, -4*y;
    0, 4*y - 1, 4*x + 4*y - 3, -4*x, 4*x, 4 - 8*y - 4*x]';
clear sqrt15
Nq = length(xhat); % numero nodi di quadratura

% Definiamo matrici con i valori che andremo ad utilizzare nell'integrale
% CASO P2
phi_matrix = zeros(Nq,Nv);
dphix_matrix = zeros(Nq,Nv);
dphiy_matrix = zeros(Nq,Nv);
for q=1:Nq
    Jphi_temp = Jphi(xhat(q), yhat(q));
    phi_matrix(q,:) = phi(xhat(q), yhat(q));
    dphix_matrix(q,:) = Jphi_temp(:,1)';
    dphiy_matrix(q,:) = Jphi_temp(:,2)';
end

for e=1:Nele
    p1 = XY(ele(e,1),:);
    p2 = XY(ele(e,2),:);
    p3 = XY(ele(e,3),:);
    area_e = geom.support.TInfo(e).Area;
    B = [(p1-p3); (p2-p3)]';
    invB = (1/det(B))*[B(2,2), -B(1,2); -B(2,1), B(1,1)];
    prodinvBinvbt = invB*invB'; % per accelerare lo calcolo una sola volta
    c = p3';
    Fe = @(x,y) c + B*[x,y]';
    for j=1:Nv
        jj = pivot(ele(e,j));
        if jj > 0
            for k=1:Nv
                kk = pivot(ele(e,k));
                if kk > 0 % assemblaggio classico di A
                    phik = phi_matrix(:, k);
                    phij = phi_matrix(:, j);
                    dphik = [dphix_matrix(:,k), dphiy_matrix(:,k)];
                    dphij = [dphix_matrix(:,j), dphiy_matrix(:,j)];
                    Djk = 0;
                    Cjk = 0;
                    Rjk = 0;
                    for q=1:Nq
                        Djk = Djk + 2*area_e*omega(q)*mu(Fe(xhat(q),yhat(q)))*dphik(q,:)*prodinvBinvbt*dphij(q,:)';
                        Cjk = Cjk + 2*area_e*omega(q)*beta(Fe(xhat(q),yhat(q)))*invB'*dphik(q,:)'*phij(q);
                        Rjk = Rjk + 2*area_e*omega(q)*sigma(Fe(xhat(q),yhat(q)))*phik(q)*phij(q);
                    end
                    A(jj,kk) = A(jj,kk) + Djk + Cjk + Rjk;
                elseif kk < 0 % assemblaggio matrice Ad legata a Dirichlet
                    phik = phi_matrix(:, k);
                    phij = phi_matrix(:, j);
                    dphik = [dphix_matrix(:,k), dphiy_matrix(:,k)];
                    dphij = [dphix_matrix(:,j), dphiy_matrix(:,j)];
                    Djk = 0;
                    Cjk = 0;
                    Rjk = 0;
                    for q=1:Nq
                        coordFe = Fe(xhat(q),yhat(q));
                        Djk = Djk + 2*area_e*omega(q)*mu(coordFe(1), coordFe(2))*dphik(q,:)*prodinvBinvbt*dphij(q,:)';
                        Cjk = Cjk + 2*area_e*omega(q)*beta(coordFe(1), coordFe(2))*invB'*dphik(q,:)'*phij(q);
                        Rjk = Rjk + 2*area_e*omega(q)*sigma(coordFe(1), coordFe(2))*phik(q)*phij(q);
                    end
                    Ad(jj,-kk) = Ad(jj,-kk) + Djk + Cjk + Rjk;
                end
            end
            for q=1:Nq
                coordFe = Fe(xhat(q),yhat(q));
                b(jj) = b(jj) + 2*area_e*omega(q)*f(coordFe(1), coordFe(2))*phij(q);
            end
        end
    end
end
% creaimo una matrice che sulle righe abbia le coordinate dei punto con
% condizioni di Dirichlet:
indexDi = geom.pivot.Di;
pointDi = XY(indexDi(:,1),:);
ud = zeros(NDi,1);
for i=1:NDi
    ud(i) = gDi(pointDi(i,1), pointDi(i,2));
end

% Imponiamo condizioni di Neuman
Ne = geom.pivot.Ne(:,1); % indice dei lati al bordo con condizioni di Ne
edgeBorders = geom.elements.borders(Ne,:,:,:);
nedgeBorders = length(Ne);
bNeumann = zeros(length(b), 1);
w = calcoloNeP2nodi();
t = [0, 0.5, 1];
for e=1:nedgeBorders
    edge = Ne(e);
    triangleBorder = edgeBorders(e,3:4);
    whereisTri = (triangleBorder > 0);
    tOB = triangleBorder(whereisTri);
    area = geom.support.TInfo(tOB).Area;
    p1 = XY(ele(tOB,1),:);
    p2 = XY(ele(tOB,2),:);
    p3 = XY(ele(tOB,3),:);
    B = [(p1-p3); (p2-p3)]';
    invB = (1/det(B))*[B(2,2), -B(1,2); -B(2,1), B(1,1)];
    c = p3';
    Fe = @(x,y) c + B*[x,y]'; % triangolo locale
    fe = @(x,y) invB*([x,y]' - c); % funzion inversa di Fe
    indexVertexTriangle = ele(tOB,:);
    indexB = geom.elements.borders(edge,1);
    indexE = geom.elements.borders(edge,2);
    indexM = geom.elements.borders(edge,5);
    indexGlobalNode = [indexB, indexM, indexE];
    % vertici globali
    Vb = XY(indexB,:);
    Ve = XY(indexE,:);
    % vertici locali
    Vbhat = fe(Vb(1), Vb(2));
    Vehat = fe(Ve(1), Ve(2));
    gamma = @(t) Vbhat + t*(Vehat - Vbhat); % trasformazione lato locale di interesse
    edgeLen = norm(Ve - Vb,2);
    for i=indexGlobalNode
        ii = pivot(i);
        if ii > 0
            integraleNe = 0;
            phiNe = @(x,y) phi(x,y)'*((indexVertexTriangle == indexB) + 0)';
            for q=1:3
                localPoint = gamma(t(q));
                globalPoint = Fe(localPoint(1), localPoint(2));
                integraleNe = integraleNe + area*edgeLen*w(q)*gNe(globalPoint(1), globalPoint(2))*phiNe(localPoint(1), localPoint(2));
            end
            b(ii) = b(ii) + integraleNe;
        end
    end
end
b = b + bNeumann;

% Risolviamo Sistema Lineare
x = A\(b-Ad*ud);
% Produciamo soluzione per output
uh = zeros(Np,1);
for j=1:Np
    jj = pivot(j);
    if jj > 0
        uh(j) = x(jj)*(jj>0) + 0;
    elseif jj < 0
        uh(j) = gDi(XY(j,1), XY(j,2));
    end
end
end