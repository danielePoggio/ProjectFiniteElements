function uh = SUPG(geom, mu, beta, f, gDi, gNe, Pk)
%SUPG for mu costante
pivot = geom.pivot.pivot;
ele = geom.elements.triangles;
XY = geom.elements.coordinates;
Np = length(XY);
Ndof = max(pivot);
NDi = -min(pivot);
Nele = length(ele);
A = zeros(Ndof,Ndof);
P = zeros(Ndof,Ndof);
G = zeros(Ndof,Ndof);
M = zeros(Ndof,1);
Ad = zeros(Ndof, NDi);
b = zeros(Ndof,1);

% definiamo i nodi di quadratura e i pesi associati
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\nodes_weights.m")
clear sqrt15
Nq = length(xhat); % numero nodi di quadratura

% Definiamo matrici con i valori che andremo ad utilizzare nell'integrale CASO P1
N1 = @(x,y) x;
N2 = @(x,y) y;
N3 = @(x,y) 1 - x - y;
if Pk == 1
    Nv = 3;
    phi = @(x,y) [N1(x,y), N2(x,y), N3(x,y)];
    Jphi = @(x,y) [1, 0; 0, 1; -1, -1];
    Hphi = zeros(2,2,6);
    mk = 1/3;
elseif Pk == 2
    Nv = 6;
    phi1 = @(x,y) 2*N1(x,y)*(N1(x,y) - 0.5);
    phi2 = @(x,y) 2*N2(x,y)*(N2(x,y) - 0.5);
    phi3 = @(x,y) 2*N3(x,y)*(N3(x,y) - 0.5);
    phi4 = @(x,y) 4*N3(x,y)*N1(x,y);
    phi5 = @(x,y) 4*N1(x,y)*N2(x,y);
    phi6 = @(x,y) 4*N2(x,y)*N3(x,y);
    phi = @(x,y) [phi1(x,y), phi2(x,y), phi3(x,y), phi4(x,y), phi5(x,y), phi6(x,y)]';
    Jphi = @(x,y) [4*x - 1, 0, 4*x + 4*y - 3, 4 - 4*y - 8*x, 4*y, -4*y;
        0, 4*y - 1, 4*x + 4*y - 3, -4*x, 4*x, 4 - 8*y - 4*x]';
    Hphi = zeros(2,2,Nv);
    Hphi(:,:,1) = [4, 0; 0,0];
    Hphi(:,:,2) = [0, 0; 0,4];
    Hphi(:,:,3) = [4, 4; 4, 4];
    Hphi(:,:,4) = [-8, -4; -4, 0];
    Hphi(:,:,5) = [0, 4; 4, 0];
    Hphi(:,:,6) = [0, -4; -4, -8];
    clear N1 N2 N3 phi1 phi2 phi3 phi4 phi5 phi6
    mk = 1/24;
end
phi_matrix = zeros(Nq,Nv);
dphix_matrix = zeros(Nq,Nv);
dphiy_matrix = zeros(Nq,Nv);
for q=1:Nq
    Jphi_temp = Jphi(xhat(q), yhat(q));
    phi_matrix(q,:) = phi(xhat(q), yhat(q));
    dphix_matrix(q,:) = Jphi_temp(:,1)';
    dphiy_matrix(q,:) = Jphi_temp(:,2)';
end
% calcoliamo viscosità artificiale
tau = zeros(Nele,1);
for e=1:Nele
    p1 = XY(ele(e,1),:);
    p2 = XY(ele(e,2),:);
    p3 = XY(ele(e,3),:);
    area_e = geom.support.TInfo(e).Area;
    B = [(p1-p3); (p2-p3)]';
    invB = (1/det(B))*[B(2,2), -B(1,2); -B(2,1), B(1,1)];
    prodinvBinvBt = invB*invB'; % per accelerare lo calcolo una sola volta
    c = p3';
    Fe = @(x,y) c + B*[x,y]';
    % calcoliamo norma di beta su E e viscosità artificiale tau(e)
    normBetaQuad = 0;
    for q=1:Nq
        coordFe = Fe(xhat(q),yhat(q));
        normBetaQuad = normBetaQuad + 2*omega(q)*beta(coordFe(1), coordFe(2))*beta(coordFe(1), coordFe(2))';
    end
    normaBeta = sqrt(normBetaQuad);
    Pe = mk*normaBeta*sqrt(area_e)*0.5/mu(coordFe(1), coordFe(2));
    if Pe < 1
        tau(e) = mk*area_e*0.25*mu(coordFe(1), coordFe(2));
    elseif Pe >= 1
        tau(e) = sqrt(area_e)*0.5/normaBeta;
    end
    for j=1:Nv
        jj = pivot(ele(e,j));
        if jj > 0
            for k=1:Nv
                kk = pivot(ele(e,k));
                if kk > 0 % assemblaggio classico di A
                    phij = phi_matrix(:, j);
                    dphik = [dphix_matrix(:,k), dphiy_matrix(:,k)];
                    dphij = [dphix_matrix(:,j), dphiy_matrix(:,j)];
                    Hphik = Hphi(:,:,k);
%                     HphikHat =  invB'*(Hphik*invB);
%                     d2phik = HphikHat(1,1) + HphikHat(2,2);
                    d2xphik = invB(1,1)*invB(1,1)*Hphik(1,1) + invB(1,1)*invB(2,1)*(Hphik(1,2) + Hphik(2,1)) + invB(2,1)*invB(2,1)*Hphik(2,2);
                    d2yphik = invB(1,2)*invB(1,2)*Hphik(1,1) + invB(1,2)*invB(2,2)*(Hphik(1,2) + Hphik(2,1)) + invB(2,2)*invB(2,2)*Hphik(2,2);
                    d2phik = d2xphik + d2yphik;
                    Djk = 0;
                    Gjk = 0;
                    Cjk = 0;
                    Pjk = 0;
                    for q=1:Nq
                        coordFe = Fe(xhat(q),yhat(q));
                        betaq = beta(coordFe(1), coordFe(2));
                        prodbetatinvBt = betaq*invB';
                        Djk = Djk + 2*area_e*omega(q)*mu(coordFe(1), coordFe(2))*dphik(q,:)*prodinvBinvBt*dphij(q,:)';
                        Cjk = Cjk + 2*area_e*omega(q)*beta(coordFe(1), coordFe(2))*invB'*dphik(q,:)'*phij(q);
                        Gjk = Gjk + tau(e)*2*area_e*omega(q)*(prodbetatinvBt*dphik(q,:)')*(prodbetatinvBt*dphij(q,:)');
                        Pjk = Pjk + tau(e)*2*area_e*omega(q)*mu(coordFe(1), coordFe(2))*d2phik*beta(coordFe(1), coordFe(2))*(invB'*dphij(q,:)');
                    end
                    A(jj,kk) = A(jj,kk) + Djk + Cjk + Gjk + Pjk;
                elseif kk < 0 % assemblaggio matrice Ad legata a Dirichlet
                    phij = phi_matrix(:, j);
                    dphik = [dphix_matrix(:,k), dphiy_matrix(:,k)];
                    dphij = [dphix_matrix(:,j), dphiy_matrix(:,j)];
                    Hphik = invB'*Hphi(:,:,k)*invB;
                    d2phik = Hphik(1,1) + Hphik(2,2);
                    Djk = 0;
                    Gjk = 0;
                    Cjk = 0;
                    Pjk = 0;
                    for q=1:Nq
                        coordFe = Fe(xhat(q),yhat(q));
                        betaq = beta(coordFe(1), coordFe(2));
                        prodbetatinvBt = betaq*invB';
                        Djk = Djk + 2*area_e*omega(q)*mu(coordFe(1), coordFe(2))*dphik(q,:)*prodinvBinvBt*dphij(q,:)';
                        Cjk = Cjk + 2*area_e*omega(q)*beta(coordFe(1), coordFe(2))*invB'*dphik(q,:)'*phij(q);
                        Gjk = Gjk + tau(e)*2*area_e*omega(q)*(prodbetatinvBt*dphik(q,:)')*(prodbetatinvBt*dphij(q,:)');
                        Pjk = Pjk + tau(e)*2*area_e*omega(q)*mu(coordFe(1), coordFe(2))*d2phik*beta(coordFe(1), coordFe(2))*(invB'*dphij(q,:)');
                    end
                    Ad(jj,-kk) = Ad(jj,-kk) + Djk + Cjk + Gjk + Pjk;
                end
            end
            for q=1:Nq
                coordFe = Fe(xhat(q),yhat(q));
                b(jj) = b(jj) + 2*area_e*omega(q)*f(coordFe(1), coordFe(2))*phij(q);
                M(jj) = M(jj) + tau(e)*2*area_e*omega(q)*f(coordFe(1), coordFe(2))*beta(coordFe(1), coordFe(2))*invB'*dphij(q,:)';
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
bNeumann = zeros(length(b), 1);
for e=1:nedgeBorders
    edge = Ne(e);
    indexB = geom.elements.borders(edge,1);
    indexE = geom.elements.borders(edge,2);
    Vb = XY(indexB,:);
    Ve = XY(indexE,:);
    edgeLen = norm(Ve - Vb,2);
    ii = geom.pivot.pivot(indexB);
    if ii > 0
        bNeumann(ii) = bNeumann(ii) + (gNe(Vb(1),Vb(2))/3 + gNe(Ve(1),Ve(2))/6)*edgeLen;
    end
    ii = geom.pivot.pivot(indexE);
    if ii > 0
        bNeumann(ii) = bNeumann(ii) + (gNe(Vb(1),Vb(2))/6 + gNe(Ve(1),Ve(2))/3)*edgeLen;
    end
end
b = b + bNeumann + M;

% Risolviamo Sistema Lineare
x = (A+P+G)\(b-Ad*ud);
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