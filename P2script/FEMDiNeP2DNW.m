function uh = FEMDiNeP2DNW(geom, mu, beta, sigma, f, gDi, gNe)
%Assemblaggio FEM Dirichlet non omogeneo
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
run("functionP2")
clear sqrt15
Nq = length(xhat); % numero nodi di quadratura

% Definiamo matrici con i valori che andremo ad utilizzare nell'integrale
% CASO P2
phi_matrix = zeros(Nq,Nv);
dphix_matrix = zeros(Nq,Nv);
dphiy_matrix = zeros(Nq,Nv);
gradphiTensor = zeros(Nv,Nq,2,1);

for q=1:Nq
    Jphi_temp = Jphi(xhat(q), yhat(q));
    phi_matrix(q,:) = phi(xhat(q), yhat(q));
    dphix_matrix(q,:) = Jphi_temp(:,1)';
    dphiy_matrix(q,:) = Jphi_temp(:,2)';
    gradphiTensor(:,q,:) = reshape(Jphi_temp, Nv,1,2);

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
    a = zeros(Nq,1);
    d = zeros(Nq,2);
    c = zeros(Nq,1);
    fq = zeros(Nq,1);
    for q=1:Nq
        coordFeq = Fe(xhat(q),yhat(q));
        a(q) = mu(coordFeq(1), coordFeq(2))*omega(q);
        d(q,:) = beta(coordFeq(1), coordFeq(2))'*omega(q);
        c(q) = sigma(coordFeq(1), coordFeq(2))*omega(q);
        fq(q) = f(coordFeq(1), coordFeq(2))*omega(q);
    end
    for j=1:Nv
        jj = pivot(ele(e,j));
        if jj > 0
            for k=1:Nv
                kk = pivot(ele(e,k));
                if kk > 0 % assemblaggio classico di A
                    gradphij = reshape(gradphiTensor(j,:,:), Nq ,2 ,1);
                    gradphik = reshape(gradphiTensor(k,:,:), Nq ,2 ,1);
                    bwp = d.*reshape(phi_matrix(:,k),Nq,1);
                    Djk = 2*area_e*a'*(gradphik.*(gradphij*prodinvBinvbt')*ones(2,1));
                    Cjk = 2*area_e*sum(bwp.*(gradphik*invB),"all");
                    Rjk = 2*area_e*sum(c.*phi_matrix(:,k).*phi_matrix(:,j), "all");
                    A(jj,kk) = A(jj,kk) + Djk + Cjk + Rjk;
                elseif kk < 0 % assemblaggio matrice Ad legata a Dirichlet
                    gradphij = reshape(gradphiTensor(j,:,:), Nq ,2 ,1);
                    gradphik = reshape(gradphiTensor(k,:,:), Nq ,2 ,1);
                    bwp = d.*reshape(phi_matrix(:,k),Nq,1);
                    Djk = 2*area_e*a'*(gradphik.*(gradphij*prodinvBinvbt')*ones(2,1));
                    Cjk = 2*area_e*sum(bwp.*(gradphik*invB),"all");
                    Rjk = 2*area_e*sum(c.*phi_matrix(:,k).*phi_matrix(:,j), "all");
                    Ad(jj,-kk) = Ad(jj,-kk) + Djk + Cjk + Rjk;
                end
            end
            b = b + sum(fq,"all");
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