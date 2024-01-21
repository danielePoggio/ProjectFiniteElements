clear all
close all
clc
 %% Definisce mesh sul dominio e controllo coincidenza del lati sul bordo
area = 0.01;
omega1Vertex = [0 1;
    1 1;
    1 2;
    0 2];
geomOmega1 = TriangolatorDomainVertex(area, omega1Vertex);
omega2Vertex = [0 0;
    2 0;
    2 1;
    0,1];
geomOmega2 = TriangolatorDomainVertex(area, omega2Vertex);

XY1 = geomOmega1.elements.coordinates;
index1 = (XY1(:,1) <= 1) & (XY1(:,2) == 1);
pointOnEdge1 = XY1(index1,:);

XY2 = geomOmega2.elements.coordinates;
index2 = (XY2(:,1) <= 1) & (XY2(:,2) == 1);
pointOnEdge2 = XY2(index2,:);

disp(sort(pointOnEdge2(:,1)) == sort(pointOnEdge1(:,1)))

N

%% Condizioni di Dirichlet nulle
% definizione del problema e della soluzione
u = @(x,y) 16*x*(1-x)*y*(1-y);
run("C:\Users\39334\Desktop\Poli\Metodi Numerici PDE\LAIB\ProjectFiniteElements\calculateDerivate.m")
gradu = @(x,y) gradu(x,y)';
d2u = @(x,y) [1,0]*Hu(x,y)*[1,0]'+ [0,1]*Hu(x,y)*[0,1]';
mu = @(x,y) 1;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(x,y) -mu(x,y)*d2u(x,y)+beta(x,y)*gradu(x,y)+sigma(x,y)*u(x,y);
n = [0,-1]'; % direzione uscente da lato su y = 0
gNe = @(x,y) mu(x,y)*n'*gradu(x,0);
gDi = @(x,y) u(x,y);

% scrivo le funzioni della base rispetto al lato (diventano funzioni
% 1D)
t = linspace(0,1,3);
phiL1 = @(t) (1-t);
phiL2 = @(t) t;
phiL = @(t) [phiL1(t), phiL2(t)]';
% sfrutto questa forma matriciale per andare a calcolare i vari
% integrali : int(phii*phij)
phiM = @(t) phiL(t)*phiL(t)';
phiTensor = zeros(2,2,3);
w = nodiQuadratura1D(3);
for k=1:3
    phiTensor(:,:,k) = w(k)*phiM(t(k));
end
phiMatrix = sum(phiTensor,3);