function checkFlag = checkOrderVertex(elementsTable, XY)
%checkOrderVertex controllo ordine vertici in triangoli
% Detailed explanation goes here
sizeT = size(elementsTable);
Nele = sizeT(1);
Nv = sizeT(2);
checkFlag = ones(Nv-3,1); % True
for e=1:Nele
    % facciamo check per ciascun triangolo 
    if sum(XY(elementsTable(e,4),:)== 0.5*(XY(elementsTable(e,1),:)+XY(elementsTable(e,2),:))) ~= 2
        checkFlag(1) = 0;
    end
    if sum(XY(elementsTable(e,5),:) == 0.5*(XY(elementsTable(e,2),:)+XY(elementsTable(e,3),:))) ~= 2
        checkFlag(2) = 0;
    end
    if sum(XY(elementsTable(e,6),:) == 0.5*(XY(elementsTable(e,1),:)+XY(elementsTable(e,3),:))) ~= 2
        checkFlag(3) = 0;
    end
end
end         