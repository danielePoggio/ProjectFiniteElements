function [uh, condA] = FEMDiNeQuadratura(geom, mu, beta, sigma, f, gDi, gNe)
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
bf = zeros(Ndof,1);

% definiamo i nodi di quadratura e i pesi associati
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\nodes_weights.m")
clear sqrt15
Nq = length(xhat); % numero nodi di quadratura

% Definiamo matrici con i valori che andremo ad utilizzare nell'integrale CASO P1
phi = @(x,y) [x, y, 1-x-y];
Jphi = @(x,y) [1, 0; 0, 1; -1, -1];
phi_matrix = zeros(Nq,3);
dphix_matrix = zeros(Nq,3);
dphiy_matrix = zeros(Nq,3);
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
    for j=1:3
        jj = pivot(ele(e,j));
        if jj > 0
            for k=1:3
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
                        coordFe = Fe(xhat(q),yhat(q));
                        Djk = Djk + 2*area_e*omega(q)*mu(coordFe(1), coordFe(2))*dphik(q,:)*prodinvBinvbt*dphij(q,:)';
                        Cjk = Cjk + 2*area_e*omega(q)*beta(coordFe(1), coordFe(2))*invB'*dphik(q,:)'*phij(q);
                        Rjk = Rjk + 2*area_e*omega(q)*sigma(coordFe(1), coordFe(2))*phik(q)*phij(q);
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
                bf(jj) = bf(jj) + 2*area_e*omega(q)*f(coordFe(1), coordFe(2))*phij(q);
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
nedgeBorders = length(edgeBorders);
bNe = zeros(Ndof, 1);
for e=1:nedgeBorders
    edge = Ne(e);
    indexB = geom.elements.borders(edge,1);
    indexE = geom.elements.borders(edge,2);
    Vb = XY(indexB,:);
    Ve = XY(indexE,:);
    edgeLen = norm(Ve - Vb,2);
    ii = geom.pivot.pivot(indexB);
    if ii > 0
        bNe(ii) = bNe(ii) + (gNe(Vb(1),Vb(2))/3 + gNe(Ve(1),Ve(2))/6)*edgeLen;
    end
    ii = geom.pivot.pivot(indexE);
    if ii > 0
        bNe(ii) = bNe(ii) + (gNe(Vb(1),Vb(2))/6 + gNe(Ve(1),Ve(2))/3)*edgeLen;
    end
end
b = bf + bNe;

condA = cond(A);
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