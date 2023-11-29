function [outputArg1,outputArg2] = integralQuadratura(mu, beta, sigma, phik, phij, geom, indexElement, B, Binvt)
%INTEGRALQUADRATURA Summary of this function goes here
%   Detailed explanation goes here
% definiamo i nodi di quadratura e i pesi associati
run("ProjectFiniteElement\funzioni\nodes_weights.m")
clear sqrt15
Nq = length(xhat); % numero nodi di quadratura
Djk = 0;
Cjk = 0;
Rjk = 0;
area_e = geom.support.TInfo(indexElement).Area; % area
% Estraggo indici globali elemento e
indexVertexElement = geom.elements.triangles(indexElement,:); % la considero Ge(k^)
% Costruisco la funzione Fe(x^,y^) = c + B*(x^,y^)'
vertexElement = XY(indexVertexElement,:);
%B = [vertexElement(2,:) - vertexElement(3,:); vertexElement(1,:) - vertexElement(3,:)];
p1 = XY(ele(indexElement,1),:);
p2 = XY(ele(indexElement,2),:);
p3 = XY(ele(indexElement,3),:);
dx1 = p3(1) - p2(1);
dx2 = p1(1) - p3(1);
dy1 = p2(2) - p3(2);
dy2 = p3(2) - p1(2);
B = [dx2, -dx1; -dy2, dy1];
invB = inv(B);
c = vertexElement(3,:)';
Fe = @(x,y) c + B*[x,y]';
for q=1:Nq
    Djk = Djk + 2*area_e*omega(q)*mu(Fe(xhat(q),yhat(q)))*dphik*invB*invB'*dphij;
    Cjk = Cjk + 2*area_e*omega(q)*beta(Fe(xhat(q),yhat(q)))'*invB'*dphik*phij;
    Rjk = Rjk + 2*area_e*omega(q)*sigma(Fe(xhat(q),yhat(q)))*phik*phij;
end
end

