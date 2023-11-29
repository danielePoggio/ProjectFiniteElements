function uh = FEMDiNoHom(geom, mu, beta, sigma, f, gDi)
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
for e=1:Nele
    p1 = XY(ele(e,1),:);
    p2 = XY(ele(e,2),:);
    p3 = XY(ele(e,3),:);
    pbar = (1/3)*(p1+p2+p3);
    mubar = mu(pbar(1),pbar(2));
    betabar = beta(pbar(1),pbar(2));
    sigmabar = sigma(pbar(1),pbar(2));
    fbar = f(pbar(1),pbar(2));
    area = areatriangolo(p1,p2,p3);
    dx1 = p3(1) - p2(1);
    dx2 = p1(1) - p3(1);
    dx3 = p2(1) - p1(1);
    dy1 = p2(2) - p3(2);
    dy2 = p3(2) - p1(2);
    dy3 = p1(2) - p2(2);
    dx = [dx1, dx2, dx3];
    dy = [dy1, dy2, dy3];
    
    for j=1:3
        jj = pivot(ele(e,j));
        if jj > 0
            for k=1:3
                kk = pivot(ele(e,k));
                if kk > 0 % assemblaggio classico di A
                    Djk = (0.25*mubar/area)*(dy(k)*dy(j) + dx(k)*dx(j));
                    Cjk = (1/6)*(betabar(1)*dy(k) + betabar(2)*dx(j));
                    Rjk = (1/12)*sigmabar*area*(1 + (j==k));
                    A(jj,kk) = A(jj,kk) + Djk + Cjk + Rjk;
                elseif kk < 0 % assemblaggio matrice Ad legata a Dirichlet
                    Djk = (0.25*mubar/area)*(dy(k)*dy(j) + dx(k)*dx(j));
                    Cjk = (1/6)*(betabar(1)*dy(k) + betabar(2)*dx(j));
                    Rjk = (1/12)*sigmabar*area*(1 + (j==k));
                    Ad(jj,-kk) = Ad(jj,-kk) + Djk + Cjk + Rjk;
                end
            end
            b(jj) = b(jj) + fbar*area/3;
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
x = A\(b-Ad*ud);
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