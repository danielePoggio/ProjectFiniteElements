function [errorL2, errorH1] = errorFunction(geom, u, gradu, uh)
% calcolare stima errore conoscendo approssimazione uh e soluzione esatta u
XY = geom.elements.coordinates;
Nele = geom.nelements.nTriangles;
ele = geom.elements.triangles;
Ne = length(geom.elements.triangles(1,:)); % numero di per ciascun elemento

% definiamo i nodi di quadratura e i pesi associati
run("ProjectFiniteElement\funzioni\nodes_weights.m")
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

errorL2 = 0;
errorH1 = 0;
for e=1:Nele
    area_e = geom.support.TInfo(e).Area; % area
    % Estraggo indici globali elemento e
    indexVertexElement = geom.elements.triangles(e,:); % la considero Ge(k^)
    % Costruisco la funzione Fe(x^,y^) = c + B*(x^,y^)'
    vertexElement = XY(indexVertexElement,:);
    %B = [vertexElement(2,:) - vertexElement(3,:); vertexElement(1,:) - vertexElement(3,:)];
    p1 = XY(ele(e,1),:);
    p2 = XY(ele(e,2),:);
    p3 = XY(ele(e,3),:);
    dx1 = p3(1) - p2(1);
    dx2 = p1(1) - p3(1);
    dy1 = p2(2) - p3(2);
    dy2 = p3(2) - p1(2);
    B = [dx2, -dx1; -dy2, dy1];
    Binvt = inv(B)';
    c = vertexElement(3,:)';
    Fe = @(x,y) c + B*[x,y]';
    integralTriangle = 0;
    temp = 0;
    for q=1:Nq
        approx_uh = 0;
        approx_grad = zeros(2,1);
        for k=1:Ne
            approx_uh = approx_uh + uh(indexVertexElement(k))*phi_matrix(q,k);
            approx_grad = approx_grad + uh(indexVertexElement(k))*[dphix_matrix(q,k), dphiy_matrix(q,k)]';
        end
        coordFe = Fe(xhat(q),yhat(q));
        integralTriangle = integralTriangle + omega(q)*2*area_e*(u(coordFe(1), coordFe(2)) - approx_uh)^2;
    end

    errorL2 = errorL2 + integralTriangle;
    vec = (gradu(coordFe(1), coordFe(2)) - Binvt*approx_grad);
    %temp = temp + 2*omega(q)*area_e*norm(vec,2)^2;
    temp = temp + 2*omega(q)*area_e*vec'*vec;
    errorH1 = errorH1 + temp;
end
errorL2 = sqrt(errorL2);
errorH1 = sqrt(errorH1);
end