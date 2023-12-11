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

% Imponiamo condizioni di Neuman
Ne = geom.pivot.Ne(:,1); % indice dei lati al bordo con condizioni di Ne
edgeBorders = geom.elements.borders(Ne,:,:,:);
nedgeBorders = length(Ne);
bNeumannP2 = zeros(Ndof, 1);
indexLocal = 1:6;
w = calcoloNeP2nodi();
t = [0, 0.5, 1];
for e=1:nedgeBorders
    % individuiamo triangolo al quale appartiene triangolo
    edge = Ne(e);
    triangleBorder = edgeBorders(e,3:4);
    whereisTri = (triangleBorder > 0);
    Tri = triangleBorder(whereisTri);
    % ricaviamo informazioni utili per trasformazioni affini
    p1 = XY(ele(Tri,1),:);
    p2 = XY(ele(Tri,2),:);
    p3 = XY(ele(Tri,3),:);
    B = [(p1-p3); (p2-p3)]';
    invB = (1/det(B))*[B(2,2), -B(1,2); -B(2,1), B(1,1)];
    c = p3';
    Fe = @(x,y) c + B*[x,y]'; % triangolo locale
    fe = @(x,y) invB*([x,y]' - c); % funzion inversa di Fe
    % vertici globali
    indexB = geom.elements.borders(edge,1);
    indexE = geom.elements.borders(edge,2);
    indexM = geom.elements.borders(edge,5);
    indexGlobalNode = [indexB, indexM, indexE];
    Vb = XY(indexB,:);
    Ve = XY(indexE,:);
    Vm = XY(indexM,:);
    % vertici locali
    Vbhat = fe(Vb(1), Vb(2));
    Vmhat = fe(Vm(1), Vm(2));
    Vehat = fe(Ve(1), Ve(2));
    for i=1:Nv
        ii = pivot(ele(Tri,i));
        if ii > 0 % grado di libertà da calcolare
            locPhi = zeros(Nv,1);
            locPhi(i) = 1;
            IntegraleLatoNe = gNe(Vb(1),Vb(2))*(phi(Vbhat(1), Vbhat(2))'*locPhi)/6 + 2*gNe(Vm(1),Vm(2))*(phi(Vmhat(1), Vmhat(2))'*locPhi)/3 + gNe(Ve(1),Ve(2))*(phi(Vehat(1), Vehat(2))'*locPhi)/6;
            bNeumannP2(ii) = bNeumannP2(ii) + gNe(Vb(1),Vb(2))*(phi(Vbhat(1), Vbhat(2))'*locPhi)/6 + 2*gNe(Vm(1),Vm(2))*(phi(Vmhat(1), Vmhat(2))'*locPhi)/3 + gNe(Ve(1),Ve(2))*(phi(Vehat(1), Vehat(2))'*locPhi)/6;
        end
    end
end