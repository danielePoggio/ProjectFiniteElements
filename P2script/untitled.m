% Imponiamo condizioni di Neuman
Ne = geom.pivot.Ne(:,1); % indice dei lati al bordo con condizioni di Ne
edgeBorders = geom.elements.borders(Ne,:,:,:);
nedgeBorders = length(Ne);
bNeumann = zeros(length(b), 1);
indexLocal = 1:6;
w = calcoloNeP2nodi();
t = [0, 0.5, 1];
for e=1:nedgeBorders
    edge = Ne(e);
    triangleBorder = edgeBorders(e,3:4);
    whereisTri = (triangleBorder > 0);
    tOB = triangleBorder(whereisTri);
    area = geom.support.TInfo(tOB).Area;
    p1 = XY(ele(tOB,1),:);
    p2 = XY(ele(tOB,2),:);
    p3 = XY(ele(tOB,3),:);
    B = [(p1-p3); (p2-p3)]';
    invB = (1/det(B))*[B(2,2), -B(1,2); -B(2,1), B(1,1)];
    c = p3';
    Fe = @(x,y) c + B*[x,y]'; % triangolo locale
    fe = @(x,y) invB*([x,y]' - c); % funzion inversa di Fe
    indexVertexTriangle = ele(tOB,:);
    indexB = geom.elements.borders(edge,1);
    localIndexB = ((indexVertexTriangle == indexB)+0)*indexLocal';
    indexE = geom.elements.borders(edge,2);
    localIndexE = ((indexVertexTriangle == indexE)+0)*indexLocal';
    indexM = geom.elements.borders(edge,5);
    localIndexM = ((indexVertexTriangle == indexM)+0)*indexLocal';
    indexGlobalNode = [indexB, indexM, indexE];
    indexLocalNode = [localIndexB, localIndexM, localIndexE];

    % vertici globali
    Vb = XY(indexB,:);
    Ve = XY(indexE,:);
    Vm = XY(indexM,:);
    % vertici locali
    Vbhat = fe(Vb(1), Vb(2));
    Vmhat = fe(Vm(1), Vm(2));
    Vehat = fe(Ve(1), Ve(2));

    gamma = @(t) Vbhat + t*(Vehat - Vbhat); % trasformazione lato locale di interesse
    edgeLen = norm(Ve - Vb,2);
    for i=indexGlobalNode
        ii = pivot(i);
%         indexVertex = geom.elements.borders(edge,i);
%         ii = pivot(indexVertex);
        if ii > 0
            indexLocalVi = indexLocal(indexVertexTriangle == indexB);
            integraleNe = 0;
            phiNe = @(x,y) phi(x,y)'*((indexVertexTriangle == indexB) + 0)';
            for q=1:3
                localPoint = gamma(t(q));
                globalPoint = Fe(localPoint(1), localPoint(2));
                integraleNe = integraleNe + area*edgeLen*w(q)*gNe(globalPoint(1), globalPoint(2))*phiNe(localPoint(1), localPoint(2));
            end
            b(ii) = b(ii) + integraleNe;
        end
    end
end