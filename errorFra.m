function [E0, E1] = errorFra(u, gradu, uh, geom, Pk)
run("nodes_weights.m")
clear sqrt15
Nq = length(xhat);
if Pk == 1
    Nv = 3;
    Phi_hat_matrix = zeros(Nq, Nv); % le funzioni di base per i P2 sono 6
    Phi_hat_matrix(:,1) = xhat';
    Phi_hat_matrix(:,2) = yhat';
    Phi_hat_matrix(:,3) = ones(Nq,1) - xhat' - yhat';

    Dx_Phi_hat_matrix = zeros(Nq,Nv);
    Dx_Phi_hat_matrix(:,1) = ones(Nq,1);
    Dx_Phi_hat_matrix(:,2) = zeros(Nq,1);
    Dx_Phi_hat_matrix(:,3) = -ones(Nq,1);

    Dy_Phi_hat_matrix = zeros(Nq, Nv);
    Dy_Phi_hat_matrix(:,1) = zeros(length(xhat), 1);
    Dy_Phi_hat_matrix(:,2) = ones(Nq,1);
    Dy_Phi_hat_matrix(:,3) = -ones(Nq,1);

elseif Pk == 2
    Nv = 6;
    what = ones(Nq, 1) - xhat' - yhat'; % N_3
    Phi_hat_matrix = zeros(Nq, Nv); % le funzioni di base per i P2 sono 6
    Phi_hat_matrix(:,1) = 2*xhat'.*(xhat'-1/2*ones(Nq, 1));
    Phi_hat_matrix(:,2) = 2*yhat'.*(yhat'-1/2*ones(Nq, 1));
    Phi_hat_matrix(:,3) = 2*what.*(what-1/2*ones(Nq, 1));
    Phi_hat_matrix(:,4) = 4*xhat'.*what;
    Phi_hat_matrix(:,5) = 4*xhat'.*yhat';
    Phi_hat_matrix(:,6) = 4*yhat'.*what;

    Dx_Phi_hat_matrix = zeros(Nq, Nv);
    Dx_Phi_hat_matrix(:,1) = 2*(xhat'-1/2*ones(length(xhat), 1)) + 2*xhat';
    Dx_Phi_hat_matrix(:,2) = zeros(length(xhat), 1);
    Dx_Phi_hat_matrix(:,3) = -2*(what-1/2*ones(length(xhat), 1)) - 2*what;
    Dx_Phi_hat_matrix(:,4) = 4*what-4*xhat';
    Dx_Phi_hat_matrix(:,5) = 4*yhat';
    Dx_Phi_hat_matrix(:,6) = -4*yhat';

    Dy_Phi_hat_matrix = zeros(Nq, Nv);
    Dy_Phi_hat_matrix(:,1) = zeros(length(xhat), 1);
    Dy_Phi_hat_matrix(:,2) = 2*(yhat'-1/2*ones(length(xhat), 1)) + 2*yhat';
    Dy_Phi_hat_matrix(:,3) = -2*(what-1/2*ones(length(xhat), 1)) - 2*what;
    Dy_Phi_hat_matrix(:,4) = -4*xhat';
    Dy_Phi_hat_matrix(:,5) = 4*xhat';
    Dy_Phi_hat_matrix(:,6) = 4*what-4*yhat';

end

Dx_u = @(x,y) [1,0]*gradu(x,y);
Dy_u = @(x,y) [0,1]*gradu(x,y);


E0 = 0;
E1 = 0;

for e = 1:geom.nelements.nTriangles
    Dx = zeros(2,1);
    Dy = zeros(2,1);
    % richiamo gli elementi della matrice B, che interverrà nella
    % trasformazione affine F_E
    Dx(1) = geom.elements.coordinates(geom.elements.triangles(e,3),1) - geom.elements.coordinates(geom.elements.triangles(e,2),1);
    Dx(2) = geom.elements.coordinates(geom.elements.triangles(e,1),1) - geom.elements.coordinates(geom.elements.triangles(e,3),1);
    Dy(1) = geom.elements.coordinates(geom.elements.triangles(e,2),2) - geom.elements.coordinates(geom.elements.triangles(e,3),2);
    Dy(2) = geom.elements.coordinates(geom.elements.triangles(e,3),2) - geom.elements.coordinates(geom.elements.triangles(e,1),2);
    
    Area = geom.support.TInfo(e).Area;
    a_3 = geom.elements.coordinates(geom.elements.triangles(e,3),:)';
    B = [Dx(2), -Dx(1);
        -Dy(2), Dy(1)];
    B_mt = (1/(2*Area))*[Dy(1), Dy(2); Dx(1), Dx(2)];
    
    %%%%%
    % soluzione ristretta al triangolo --> vettore contenente i valori (della soluzione approssimata) nei 3
    % vertici del triangolo. A differenza dei P1, qui ho 6 gdl--> 6 valori
    % nei nodi
    U_T = zeros(Nv,1);
    for i = 1:Nv
        U_T(i) = uh(geom.elements.triangles(e,i));
    end
    %%%%%

    for q = 1:length(xhat)
        F_T = a_3 + B*[xhat(q), yhat(q)]'; % trasformazione affine
        E0 = E0 + omega(q)*2*Area*(u(F_T(1), F_T(2))-Phi_hat_matrix(q,:)*U_T)^2;
        E1 = E1 + omega(q)*2*Area*( [Dx_u(F_T(1), F_T(2)); Dy_u(F_T(1), F_T(2))]-B_mt*[Dx_Phi_hat_matrix(q,:)*U_T; Dy_Phi_hat_matrix(q,:)*U_T] )'*( [Dx_u(F_T(1), F_T(2)); Dy_u(F_T(1), F_T(2))]-B_mt*[Dx_Phi_hat_matrix(q,:)*U_T; Dy_Phi_hat_matrix(q,:)*U_T] );
    end
end
E0 = sqrt(E0);
E1 = sqrt(E1);
end