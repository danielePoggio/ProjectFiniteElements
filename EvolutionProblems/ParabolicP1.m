function uh = ParabolicP1(geom, deltat, Nt, rho, mu, beta, sigma, f, gDi, gNe, dtgDi, u0)
%Assemblaggio FEM P1 -> supponiamo che nodi possono variare nello spazio ma
% %non nel tempo
pivot = geom.pivot.pivot;
indexDof = (pivot > 0);
ele = geom.elements.triangles;
XY = geom.elements.coordinates;
Np = length(XY);
Ndof = max(pivot);
NDi = -min(pivot);
Nele = length(ele);
A = zeros(Ndof,Ndof);
B = zeros(Ndof,Ndof);
Ad = zeros(Ndof, NDi);
Bd = zeros(Ndof, NDi);
uh = zeros(Np, Nt+1);
for i=1:Np
    V = XY(i,:);
    uh(i,1) = u0(V(1),V(2));
end

% definiamo i nodi di quadratura e i pesi associati
run("nodes_weights.m")
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
    BFe = [(p1-p3); (p2-p3)]';
    invB = (1/det(BFe))*[BFe(2,2), -BFe(1,2); -BFe(2,1), BFe(1,1)];
    prodinvBinvbt = invB*invB'; % per accelerare lo calcolo una sola volta
    c = p3';
    Fe = @(x,y) c + BFe*[x,y]';
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
                    Bjk = 0;
                    for q=1:Nq
                        Djk = Djk + 2*area_e*omega(q)*mu(Fe(xhat(q),yhat(q)))*dphik(q,:)*prodinvBinvbt*dphij(q,:)';
                        Cjk = Cjk + 2*area_e*omega(q)*beta(Fe(xhat(q),yhat(q)))*invB'*dphik(q,:)'*phij(q);
                        Rjk = Rjk + 2*area_e*omega(q)*sigma(Fe(xhat(q),yhat(q)))*phik(q)*phij(q);
                        Bjk = Bjk + 2*area_e*omega(q)*rho(Fe(xhat(q),yhat(q)))*phik(q)*phij(q);
                    end
                    A(jj,kk) = A(jj,kk) + Djk + Cjk + Rjk;
                    B(jj,kk) = B(jj,kk) + Bjk;
                elseif kk < 0 % assemblaggio matrice Ad legata a Dirichlet
                    phik = phi_matrix(:, k);
                    phij = phi_matrix(:, j);
                    dphik = [dphix_matrix(:,k), dphiy_matrix(:,k)];
                    dphij = [dphix_matrix(:,j), dphiy_matrix(:,j)];
                    Djk = 0;
                    Cjk = 0;
                    Rjk = 0;
                    Bjk = 0;
                    for q=1:Nq
                        Djk = Djk + 2*area_e*omega(q)*mu(Fe(xhat(q),yhat(q)))*dphik(q,:)*prodinvBinvbt*dphij(q,:)';
                        Cjk = Cjk + 2*area_e*omega(q)*beta(Fe(xhat(q),yhat(q)))*invB'*dphik(q,:)'*phij(q);
                        Rjk = Rjk + 2*area_e*omega(q)*sigma(Fe(xhat(q),yhat(q)))*phik(q)*phij(q);
                        Bjk = Bjk + 2*area_e*omega(q)*rho(Fe(xhat(q),yhat(q)))*phik(q)*phij(q);
                    end
                    Ad(jj,-kk) = Ad(jj,-kk) + Djk + Cjk + Rjk;
                    Bd(jj,-kk) = Bd(jj,-kk) + Bjk;
                end
            end
        end
    end
end
for n=1:Nt
    tstep = n*deltat;
    b1 = zeros(Ndof,1);
    b2 = zeros(Ndof,1);

    for e=1:Nele
        p1 = XY(ele(e,1),:);
        p2 = XY(ele(e,2),:);
        p3 = XY(ele(e,3),:);
        area_e = geom.support.TInfo(e).Area;
        BFe = [(p1-p3); (p2-p3)]';
        c = p3';
        Fe = @(x,y) c + BFe*[x,y]';
        % Calcolo termine noto b(tstep)
        for j=1:3
            jj = pivot(ele(e,j));
            if jj > 0
                for q=1:Nq
                    coordFe = Fe(xhat(q),yhat(q));
                    b1(jj) = b1(jj) + 2*area_e*omega(q)*f(tstep, coordFe(1), coordFe(2))*phij(q);
                    b2(jj) = b2(jj) + 2*area_e*omega(q)*f(tstep+deltat, coordFe(1), coordFe(2))*phij(q);

                end
            end
        end
    end
    % creaimo una matrice che sulle righe abbia le coordinate dei punto con
    % condizioni di Dirichlet:
    indexDi = geom.pivot.Di;
    pointDi = XY(indexDi(:,1),:);
    ud1 = zeros(NDi,1);
    dtud1 = zeros(NDi,1);
    ud2 = zeros(NDi,1);
    dtud2 = zeros(NDi,1);
    for i=1:NDi
        ud1(i) = gDi(tstep, pointDi(i,1), pointDi(i,2));
        dtud1(i) = dtgDi(tstep, pointDi(i,1), pointDi(i,2));
        ud2(i) = gDi(tstep+deltat, pointDi(i,1), pointDi(i,2));
        dtud2(i) = dtgDi(tstep+deltat, pointDi(i,1), pointDi(i,2));
    end

    % Imponiamo condizioni di Neuman
    Ne = geom.pivot.Ne(:,1); % indice dei lati al bordo con condizioni di Ne
    edgeBorders = geom.elements.borders(Ne,:,:,:);
    nedgeBorders = length(edgeBorders);
    bNeumann1 = zeros(Ndof, 1);
    bNeumann2 = zeros(Ndof, 1);
    for e=1:nedgeBorders
        edge = Ne(e);
        indexB = geom.elements.borders(edge,1);
        indexE = geom.elements.borders(edge,2);
        Vb = XY(indexB,:);
        Ve = XY(indexE,:);
        edgeLen = norm(Ve - Vb,2);
        ii = geom.pivot.pivot(indexB);
        if ii > 0
            bNeumann1(ii) = bNeumann1(ii) + (gNe(tstep, Vb(1),Vb(2))/3 + gNe(tstep, Ve(1),Ve(2))/6)*edgeLen;
            bNeumann2(ii) = bNeumann2(ii) + (gNe(tstep+deltat, Vb(1),Vb(2))/3 + gNe(tstep+deltat,Ve(1),Ve(2))/6)*edgeLen;
        end
        ii = geom.pivot.pivot(indexE);
        if ii > 0
            bNeumann1(ii) = bNeumann1(ii) + (gNe(tstep, Vb(1),Vb(2))/6 + gNe(tstep, Ve(1),Ve(2))/3)*edgeLen;
            bNeumann2(ii) = bNeumann2(ii) + (gNe(tstep+deltat, Vb(1),Vb(2))/6 + gNe(tstep+deltat, Ve(1),Ve(2))/3)*edgeLen;
        end
    end
    matrix = (B+0.5*deltat*A);
    termineNoto = (B-0.5*deltat*A)*uh(indexDof,n) - 0.5*deltat*(Bd*dtud1 + Ad*ud1 + Bd*dtud2 + Ad*ud2 - (b1+b2) - (bNeumann1+bNeumann2));
    % Risolviamo Sistema Lineare
    x = matrix\termineNoto;
    % Produciamo soluzione per output
    for j=1:Np
        jj = pivot(j);
        if jj > 0
            uh(j, n+1) = x(jj)*(jj>0) + 0;
        elseif jj < 0
            uh(j,n+1) = gDi(tstep, XY(j,1), XY(j,2));
        end
    end
end
end