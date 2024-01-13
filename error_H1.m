function E1 = error_H1(gradu, uh, xhat, yhat, omega, geom)
Nq = length(xhat);
if Pk == 1
    Nv = 3;
    Dx_Phi_hat_matrix = zeros(length(xhat),Nv);
    Dx_Phi_hat_matrix(:,1) = ones(Nq,1);
    Dx_Phi_hat_matrix(:,2) = zeros(Nq,1);
    Dx_Phi_hat_matrix(:,3) = -ones(Nq,1);

    Dy_Phi_hat_matrix = zeros(length(xhat),6);
    Dy_Phi_hat_matrix(:,1) = zeros(length(xhat), 1);
    Dy_Phi_hat_matrix(:,2) = ones(Nq,1);
    Dy_Phi_hat_matrix(:,3) = -ones(Nq,1);

elseif Pk == 2
    Nv = 6;
    what = ones(length(xhat), 1) - xhat' - yhat'; % N_3
    Dx_Phi_hat_matrix = zeros(length(xhat),6);
    Dx_Phi_hat_matrix(:,1) = 2*(xhat'-1/2*ones(length(xhat), 1)) + 2*xhat';
    Dx_Phi_hat_matrix(:,2) = zeros(length(xhat), 1);
    Dx_Phi_hat_matrix(:,3) = -2*(what-1/2*ones(length(xhat), 1)) - 2*what;
    Dx_Phi_hat_matrix(:,4) = 4*what-4*xhat';
    Dx_Phi_hat_matrix(:,5) = 4*yhat';
    Dx_Phi_hat_matrix(:,6) = -4*yhat';

    Dy_Phi_hat_matrix = zeros(length(xhat),6);
    Dy_Phi_hat_matrix(:,1) = zeros(length(xhat), 1);
    Dy_Phi_hat_matrix(:,2) = 2*(yhat'-1/2*ones(length(xhat), 1)) + 2*yhat';
    Dy_Phi_hat_matrix(:,3) = -2*(what-1/2*ones(length(xhat), 1)) - 2*what;
    Dy_Phi_hat_matrix(:,4) = -4*xhat';
    Dy_Phi_hat_matrix(:,5) = 4*xhat';
    Dy_Phi_hat_matrix(:,6) = 4*what-4*yhat';

end

Dx_u = @(x,y) [1,0]*gradu(x,y);
Dy_u = @(x,y) [0,1]*gradu(x,y);

E1 = 0;

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
    B_mt = (1/(2*Area))*[Dy(1), Dy(2); Dx(1), Dx(2)];
    
    %%%%%
    % soluzione ristretta al triangolo --> vettore contenente i valori (della soluzione approssimata) nei 3
    % vertici del triangolo
    U_T = zeros(Nv,1);
    for i = 1:Nv
        U_T(i) = uh(geom.elements.triangles(e,i));
    end

    for q = 1:Nq
        F_T1 = a_3 + B*[xhat(q), yhat(q)]'; % trasformazione affine
        E1 = E1 + omega(q)*2*Area*( [Dx_u(F_T1(1), F_T1(2)); Dy_u(F_T1(1), F_T1(2))]-B_mt*[Dx_Phi_hat_matrix(q,:)*U_T; Dy_Phi_hat_matrix(q,:)*U_T] )'*( [Dx_u(F_T1(1), F_T1(2)); Dy_u(F_T1(1), F_T1(2))]-B_mt*[Dx_Phi_hat_matrix(q,:)*U_T; Dy_Phi_hat_matrix(q,:)*U_T] );
    end
end
end