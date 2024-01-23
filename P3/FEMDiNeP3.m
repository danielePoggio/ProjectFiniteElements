function [uh, idxValue] = FEMDiNeP3(geom, mu, beta, sigma, f, gDi, gNe)
%Assemblaggio FEM Dirichlet non omogeneo, calcolo coeff con nodi di
%quadratura invece che approssimare con il baricentro
pivot = geom.pivot.pivot;
ele = geom.elements.triangles;
XY = geom.elements.coordinates;
Nvertex = length(XY);
Ndeg = sum((pivot ~= 0), "all");
Ndof = geom.pivot.Ndof;
NDi = max(-min(pivot));
Nele = length(ele);
A = zeros(Ndof,Ndof);
Ad = zeros(Ndof, NDi);
b = zeros(Ndof,1);

% ricaviamo nodi di quadratura
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\nodes_weights.m")
clear sqrt15
Nq = length(xhat); % numero nodi di quadratura
% estraiamo funzioni P3
Ndegree = 10;
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\P3\P3referencebasis.m");
% Definiamo matrici con i valori che andremo ad utilizzare nell'integrale
phi_matrix = zeros(Nq,Ndegree);
dphix_matrix = zeros(Nq,Ndegree);
dphiy_matrix = zeros(Nq,Ndegree);
for i = 1:Ndegree
    n = zeros(Ndegree,1);
    n(i) = 1;
    phi_i = @(x,y) n'*phiP3(x,y);
    [gradphi_i, ~] = calculateDerivate(phi_i);
    for q=1:Nq
        phi_matrix(q,i) = phi_i(xhat(q),yhat(q));
        grad_temp = gradphi_i(xhat(q),yhat(q));
        dphix_matrix(q,i) = grad_temp(1);
        dphiy_matrix(q,i) = grad_temp(2);
    end
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
    for j=1:Ndegree
        vertexIndexj = vertexIndexP3(j);
        globalIndexj = ele(e,vertexIndexj);
        indexGdlj = localIndexP3(j);
        jj = pivot(globalIndexj, indexGdlj);
        if jj > 0
            for k=1:Ndegree
                vertexIndexk = vertexIndexP3(k);
                globalIndexk = ele(e,vertexIndexk);
                indexGdlk = localIndexP3(k);
                kk = pivot(globalIndexk, indexGdlk);
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
                else % assemblaggio matrice Ad legata a Dirichlet
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
bNeumannP2 = zeros(Ndof, 1);
t = linspace(0,1,3);
% scrivo le funzioni della base rispetto al lato (diventano funzioni
% 1D)
phiL1 = @(t) 2*(t-0.5)*(t-1);
phiL2 = @(t) -4*t*(t-1);
phiL3 = @(t) 2*t*(t-0.5);
phiL = @(t) [phiL1(t), phiL2(t), phiL3(t)]';
% sfrutto questa forma matriciale per andare a calcolare i vari
% integrali : int(phii*phij)
phiM = @(t) phiL(t)*phiL(t)';
phiTensor = zeros(3,3,3);
w = nodiQuadratura1D(3);
for k=1:3
    phiTensor(:,:,k) = w(k)*phiM(t(k));
end
phiMatrix = sum(phiTensor,3);
for e=1:nedgeBorders
    edge = Ne(e);
    indexB = geom.elements.borders(edge,1);
    indexE = geom.elements.borders(edge,2);
    indexM = geom.elements.borders(edge,5);
    Vb = XY(indexB,:);
    Ve = XY(indexE,:);
    Vm = XY(indexM,:);
    edgeLen = norm(Ve - Vb,2);
    gN = [gNe(Vb(1), Vb(2)) gNe(Vm(1), Vm(2)) gNe(Ve(1), Ve(2))]';
    finalIntegrals = edgeLen*(phiMatrix*gN);
    indexNode = [indexB, indexM, indexE];
    for i=1:3
        idx = indexNode(i);
        ii = geom.pivot.pivot(idx);
        if ii > 0
            bNeumannP2(ii) = bNeumannP2(ii) + finalIntegrals(i);
        end
    end
end
b = b + bNeumannP2;

% condA = cond(A);
% Risolviamo Sistema Lineare
x = A\(b-Ad*ud);
% Produciamo soluzione per output
uh = zeros(Ndeg,1);
localCount = 0;
idxValue = zeros(Nvertex,1);
for i=1:Nvertex
    for j = 1:3
        ii = pivot(i,j);
        if ii ~= 0 % allora è un vincolo
            localCount = localCount + 1;
            if ii > 0
                uh(localCount) = x(ii)*(jj>0) + 0;
            elseif ii < 0
                uh(localCount) = gDi(XY(i,1), XY(i,2));
            end
            if j==1
                idxValue(i) = localCount;
            end
        end
    end
end 
end

