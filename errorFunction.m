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
    Nv = 3;
    phi = @(x,y) [x, y, 1-x-y];
    Jphi = @(x,y) [1, 0; 0, 1; -1, -1];
    phi_matrix = zeros(Nq,Nv);
    dphix_matrix = zeros(Nq,Nv);
    dphiy_matrix = zeros(Nq,Nv);
    for q=1:Nq
        Jphi_temp = Jphi(xhat(q), yhat(q));
        phi_matrix(q,:) = phi(xhat(q), yhat(q));
        dphix_matrix(q,:) = Jphi_temp(:,1)';
        dphiy_matrix(q,:) = Jphi_temp(:,2)';
    end
elseif Pk == 2
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
    phi_matrix = zeros(Nq,Nv);
    dphix_matrix = zeros(Nq,Nv);
    dphiy_matrix = zeros(Nq,Nv);
    for q=1:Nq
        Jphi_temp = Jphi(xhat(q), yhat(q));
        phi_matrix(q,:) = phi(xhat(q), yhat(q));
        dphix_matrix(q,:) = Jphi_temp(:,1)';
        dphiy_matrix(q,:) = Jphi_temp(:,2)';
    end
end
errorL2 = 0;
errorH1 = 0;
for e=1:Nele
    % Estraggo indici globali elemento e
    indexVertexElement = geom.elements.triangles(e,:); % la considero Ge(k^)
    % Costruisco la funzione Fe(x^,y^) = c + B*(x^,y^)'
    p1 = XY(ele(e,1),:);
    p2 = XY(ele(e,2),:);
    p3 = XY(ele(e,3),:);
    area_e = geom.support.TInfo(e).Area;
    B = [(p1-p3); (p2-p3)]';
    invB = (1/det(B))*[B(2,2), -B(1,2); -B(2,1), B(1,1)];
    c = p3';
    Fe = @(x,y) c + B*[x,y]';
    integralTriangle = 0;
    temp = 0;
    for q=1:Nq
        approx_uh = 0;
        approx_grad = zeros(2,1);
        for k=1:Nv
            approx_uh = approx_uh + uh(indexVertexElement(k))*phi_matrix(q,k);
            approx_grad = approx_grad + uh(indexVertexElement(k))*[dphix_matrix(q,k), dphiy_matrix(q,k)]';
        end
        coordFe = Fe(xhat(q),yhat(q));
        integralTriangle = integralTriangle + omega(q)*2*area_e*(u(coordFe(1), coordFe(2)) - approx_uh)^2;
        vec = (gradu(coordFe(1), coordFe(2)) - invB'*approx_grad);
        temp = temp + 2*omega(q)*area_e*(vec'*vec);
    end
    errorL2 = errorL2 + integralTriangle;
    %temp = temp + 2*omega(q)*area_e*norm(vec,2)^2;
    errorH1 = errorH1 + temp;
end
errorL2 = sqrt(errorL2);
errorH1 = sqrt(errorH1);
end