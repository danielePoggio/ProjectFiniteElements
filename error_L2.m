function [E0] = error_L2(u, uh, geom, Pk)

run("nodes_weights.m")
clear sqrt15
Nq = length(xhat);
if Pk == 1
    Nv = 3;
    Phi_hat_matrix = zeros(Nq, Nv); % le funzioni di base per i P2 sono 6
    Phi_hat_matrix(:,1) = xhat';
    Phi_hat_matrix(:,2) = yhat';
    Phi_hat_matrix(:,3) = ones(Nq,1) - xhat' - yhat';

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

end


E0 = 0;

for e = 1:geom.nelements.nTriangles
    
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
        F_T0 = a_3 + B*[xhat(q), yhat(q)]'; % trasformazione affine
        E0 = E0 + omega(q)*2*Area*(u(F_T0(1), F_T0(2))-Phi_hat_matrix(q,:)*U_T)^2;
    end
end

end