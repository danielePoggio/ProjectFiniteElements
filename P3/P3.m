function geom = P3(geom)
% P3 Aggiormento struttura dati

% Definiamo Baricentri e aggiorniamo XY e la struttura ele
XY = geom.elements.coordinates;
Np = length(XY);
ele = geom.elements.triangles;
[Nele, Nv] = size(ele);
ele2 = zeros(Nele, Nv+1);
ele2(:, 1:Nv) = ele;
XY2 = zeros(Np+Nele,2);
XY2(1:Np,:) = XY;
pivot = geom.pivot.pivot;

Ndof = max(pivot);
pivot2 = zeros(Np + Nele,1);
pivot2(1:Np,1) = pivot;
for ntri = 1:Nele
    XY2(Np+ntri,:) = geom.support.TInfo(ntri).CG;
    ele2(ntri,Nv+1) = Np+ntri;
    Ndof = Ndof + 1;
    pivot2(Np+ntri) = Ndof;
end
geom.elements.triangles = ele2;
geom.elements.coordinates = XY2;
geom.pivot.pivot = pivot2;

% Aggiorniamo nodelist
nodelist = geom.pivot.nodelist;
nodelist2 = zeros(Np + Nele,1);
nodelist2(1:Np,1) = nodelist;
geom.pivot.nodelist = nodelist2;

% definiamo degreeFree che indica quanti gradi di libertà ci sono per
% ciascun nodo
Di = geom.pivot.Di(:,1);
degreeFree = zeros(Np+Nele, 3);
degreeFree(1:Np,:) = ones(Np,3);
degreeFree(Di, 1) = zeros(length(Di),1);
degreeFree(Np+1:Nele+Np,1) = ones(Nele,1);
Ndof = sum(degreeFree,"all");
geom.pivot.gdl = degreeFree;
geom.pivot.Ndof = Ndof;
pivot2real = zeros(Np+Nele, 3);
countPos = 0;
countNeg = 0;
for i=1:Np+Nele
    for j=1:3
        if degreeFree(i,j) == 1
            countPos = countPos + 1;
            pivot2real(i,j) = countPos;
        elseif degreeFree(i,j) == 0 && (i<Np+1)
            countNeg = countNeg - 1;
            pivot2real(i,j) = countNeg;
        end
    end
end
geom.pivot.pivot = pivot2real;
end


