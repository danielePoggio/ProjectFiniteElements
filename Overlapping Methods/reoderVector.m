function XYnew = reoderVector(XY, indexNe, location)
%reoderVector riordino vettore nell'ordine in cui voglio
XYnew = zeros(size(XY));
XYNe = XY(indexNe,:);
XYother = XY;
XYother(indexNe,:) = [];

if location == 'd'
    XYnew(1:length(XYother),:) = XYother;
    XYnew(length(XYother)+1:length(XYother)+length(XYNe),:) = XYNe;
elseif location == 't'
    XYnew(1:length(XYNe),:) = XYNe;
    XYnew(length(XYNe)+1:length(XYother)+length(XYNe),:) = XYother;
end

end