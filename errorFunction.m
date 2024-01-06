function [errorL2, errorH1] = errorFunction(geom, u, gradu, uh, Pk)
% calcolare stima errore conoscendo approssimazione uh e soluzione esatta u
XY = geom.elements.coordinates;
Nele = geom.nelements.nTriangles;
ele = geom.elements.triangles;
Nv = length(geom.elements.triangles(1,:)); % numero di per ciascun elemento

% definiamo i nodi di quadratura e i pesi associati
run("nodes_weights.m")
clear sqrt15
Nq = length(xhat); % numero nodi di quadratura

% Definiamo matrici con i valori che andremo ad utilizzare nell'integrale CASO P1
if Pk == 1
    phi = @(x,y) [x, y, 1-x-y];
    Jphi = @(x,y) [1, 0; 0, 1; -1, -1];
    phi_matrix = zeros(Nq,Nv);
    gradphi_matrix = zeros(2, Nv, Nq);
    for q=1:Nq
        phi_matrix(q,:) = phi(xhat(q), yhat(q));
        gradphi_matrix(:,:,q) = Jphi(xhat(q), yhat(q))';
    end
elseif Pk == 2
    N1 = @(x,y) x;
    N2 = @(x,y) y;
    N3 = @(x,y) 1 - x - y;
    phi1 = @(x,y) 2*N1(x,y)*(N1(x,y) - 0.5);
    phi2 = @(x,y) 2*N2(x,y)*(N2(x,y) - 0.5);
    phi3 = @(x,y) 2*N3(x,y)*(N3(x,y) - 0.5);
    phi4 = @(x,y) 4*N3(x,y)*N1(x,y);
    phi5 = @(x,y) 4*N1(x,y)*N2(x,y);
    phi6 = @(x,y) 4*N2(x,y)*N3(x,y);
    phi = @(x,y) [phi1(x,y), phi2(x,y), phi3(x,y), phi4(x,y), phi5(x,y), phi6(x,y)];

    Jphi = @(x,y) [4*x - 1, 0, 4*x + 4*y - 3, 4 - 4*y - 8*x, 4*y, -4*y;
        0, 4*y - 1, 4*x + 4*y - 3, -4*x, 4*x, 4 - 8*y - 4*x]';
    phi_matrix = zeros(Nq,Nv);
    gradphi_matrix = zeros(2, Nv, Nq);
    for q=1:Nq
        phi_matrix(q,:) = phi(xhat(q), yhat(q));
        gradphi_matrix(:,:,q) = Jphi(xhat(q), yhat(q))';
    end
end
A = zeros(Nv,Nv);
for k=1:Nv
    for j=1:Nv
        A(k,j) = A(k,j) + omega*(phi_matrix(:,k).*phi_matrix(:,j));
    end
end
errorL2 = 0;
errorH1 = 0;
for e=1:Nele
    % Estraggo indici globali elemento e
    indexVertexElement = ele(e,:); % la considero Ge(k^)
    % Costruisco la funzione Fe(x^,y^) = c + B*(x^,y^)'
    p1 = XY(ele(e,1),:);
    p2 = XY(ele(e,2),:);
    p3 = XY(ele(e,3),:);
    area_e = geom.support.TInfo(e).Area;
    B = [(p1-p3); (p2-p3)]';
    invB = (1/det(B))*[B(2,2), -B(1,2); -B(2,1), B(1,1)];
    c = p3';
    Fe = @(x,y) c + B*[x,y]';
    % calcoliamo valori di u e gradu nei gradi di libertà del triangolo
    uLocal = zeros(Nv,1);
    graduLocal = zeros(2,Nv);
    for k=1:Nv
        coord = XY(ele(e,k),:); 
        uLocal(k) = uh(ele(e,k));
        graduLocal(:,k) = gradu(coord(1),coord(2));
    end
    % Calcoliamo errore norma L2 su elemento E
    uPhi = zeros(Nv,1);
    gradUgradPhi = zeros(Nv,1);
    gradPhigradU = zeros(Nv,1);
    gradgrad = zeros(Nv,Nv);
    u_ = zeros(Nq,1); % valore di u nei nodi di interpolazione
    gradu_ = zeros(2,Nq);
    uSquared = zeros(Nq,1);
    graduSquared = zeros(Nq,1);
    for q=1:Nq
        coordFe = Fe(xhat(q),yhat(q));
        uq = u(coordFe(1), coordFe(2));
        graduq = gradu(coordFe(1), coordFe(2));
        uSquared(q) = uq*uq;
        graduSquared(q) = graduq'*graduq;
        u_(q) = uq;
        gradu_(:, q) = graduq;
    end
    normauSquared = omega*uSquared;
    normaGraduSquared = omega*graduSquared;
    for k=1:Nv
        uPhi(k) = omega*(u_.*phi_matrix(:,k));
        gradUgradPhi(k) = sum(gradu_.*(invB'*reshape(gradphi_matrix(:,k,:), 2,7)),1)*omega' + sum(reshape(gradphi_matrix(:,k,:), 2,7).*(invB*gradu_),1)*omega';
        for j = 1:Nv
            gradgrad(k,j) = sum(reshape(gradphi_matrix(:,k,:), 2,7).*(invB*(invB'*reshape(gradphi_matrix(:,j,:), 2,7))),1)*omega';
        end
    end
    errorL2 = errorL2 + 2*area_e*(normauSquared - 2*uLocal'*uPhi + uLocal'*A*uLocal);
    errorH1 = errorH1 + 2*area_e*(normaGraduSquared - uLocal'*(gradUgradPhi) + uLocal'*gradgrad*uLocal);
end
errorL2 = sqrt(errorL2);
errorH1 = sqrt(errorH1);
end