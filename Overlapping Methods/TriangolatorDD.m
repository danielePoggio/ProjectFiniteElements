function geom = TriangolatorDD(area)
if(~exist('bbtr30'))
     addpath('../bbtr30')
     disp('../bbtr30 added to the path')
end

%----------------------------------------------------------------------------
%
% Triangolazione di un dominio quadrata
% con condizioni di Dirichlet sul bordo
%
%----------------------------------------------------------------------------
%
%  Autore: Stefano Berrone
%  Politecnico di Torino
%
%----------------------------------------------------------------------------

% clc
% clear all

% -------------------------------
% Inserimento dei vertici
% -------------------------------

Domain.InputVertex = [ 0 0
5 0
5 5
0 5
1 1
4 1.25
1 1.5
4 4
% 4 3.5
% 1 3.75
3 0
3 5];
Domain.Boundary.Values = [1 11 2 3 10 4] ;
Domain.Holes.Hole = [];
% Domain.Holes.Hole(1).Values = 5:7;
% Domain.Holes.Hole(2).Values = 8:10;
% Domain.Segments.Segment(1).Values = [4 10];
% Domain.Segments.Segment(2).Values = [10 6];
% Domain.Segments.Segment(3).Values = [6 2];
Domain.Segments.Segment(4).Values = [11 12];
BC.Values = [2 3 4 6 7 8];
BC.Boundary.Values = [1 0 2 1 0 2];
BC.Holes.Hole(1).Values = [3 3 3];
BC.Holes.Hole(2).Values = [4 4 4];
BC.Segments.Segment(1).Values = [0];
BC.Segments.Segment(2).Values = [0];

% --------------------------------------------
% Creazione della triangolazione e plottaggio
% --------------------------------------------

BC.Segments.Segment(3).Values = [0];
BC.Segments.Segment(4).Values = [0];
BC.InputVertexValues = [1 1 0 0 3 3 3 0 0 0 ];
RefiningOptions.CheckArea = 'Y';
RefiningOptions.CheckAngle = 'N';
RefiningOptions.AreaValue = area;
RefiningOptions.AngleValue = [ ];
RefiningOptions.Subregions = [ ];
[geom] = bbtr30(Domain,BC,RefiningOptions);
draw_grid (geom,1);

% --------------------------------------------------
% --------------------------------------------------


% --------------------------------------------------
% Rielaborazione dei prodotti del triangolatore
% per un piu` agevole trattamento delle condizioni
% al contorno
% --------------------------------------------------

geom.elements.coordinates = geom.elements.coordinates(...
				1:geom.nelements.nVertexes,:);
geom.elements.triangles = geom.elements.triangles(...
				1:geom.nelements.nTriangles,:);
geom.elements.borders = geom.elements.borders(...
				1:geom.nelements.nBorders,:);
geom.elements.neighbourhood = geom.elements.neighbourhood(...
				1:geom.nelements.nTriangles,:);

% --------------------------------------------------

j  = 1;
Dj = 1;
for i=1:size(geom.pivot.nodelist)
     if geom.pivot.nodelist(i)==0
        geom.pivot.pivot(i)=j;
        j = j+1;
     else
        geom.pivot.pivot(i)=-Dj;
        Dj = Dj + 1;
     end
end

% --------------------------------------------------

geom.pivot.pivot = transpose(geom.pivot.pivot);

% --------------------------------------------------

% geom.pivot.Di dopo le operazioni seguenti contiene l`indice dei nodi
% di Dirichlet e il corrispondente marker

[X,I] = sort(geom.pivot.Di(:,1));
geom.pivot.Di = geom.pivot.Di(I,:);

clear X I;
end