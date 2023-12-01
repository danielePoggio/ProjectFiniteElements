
nnode = geom.nelements.nVertexes;

for l=1:geom.nelements.nBorders %numero lati triangolazione

  n(1)=geom.elements.borders(l,1); % edge starting node
  n(2)=geom.elements.borders(l,2); % edge ending node
  e(1)=geom.elements.borders(l,3); % side element triangolo da una parte
  e(2)=geom.elements.borders(l,4); % side element triang altra parte 

  %per ogni lato aggiungo un nodo
  nnode = nnode + 1;
%trovo punto medio del lato che aggiungo alla triangolazione
  geom.elements.coordinates(nnode,:) = (geom.elements.coordinates(n(1),:)+...
					geom.elements.coordinates(n(2),:))/2; % coordinates of the edge midpoint
%associo ad ogni lato il suo punto medio
%creo una quinta colonna in cui ho le coordinate del punto medio
  geom.elements.borders(l,5)=nnode; % to connect the edge with its midpoint
  
  idx = [1 2 3];

  for el=e
    
    if(el ~= -1) %escludo il caso in cui il lato è condiviso da un solo triangolo
        %ora devo capire dove devo aggiungere il nuovo indice del nodo
      acc = 0;
      acc = idx * ( geom.elements.triangles(el,1:3)==n(1) )'; %quale dei tre vertici ha indice globale n(1) inizio lato
      %acc avrà due 1 e tutti 0
      acc = acc + idx * ( geom.elements.triangles(el,1:3)==n(2) )';

      switch acc
    case 3 %somma indici locali =3 (tra vertice 1 e 2)
	  geom.elements.triangles(el,4) = nnode;
	case 4
	  geom.elements.triangles(el,6) = nnode;
	case 5
	  geom.elements.triangles(el,5) = nnode;
	otherwise
	  disp('sconoscuto');
      end % switch acc      
    end %qui ho considerato tutti i lati ma non ho considerato le condizioni al contorno non omogenee

  end % for el=e
%da info locali voglio ricostruire la posizione sul dominio per imporre la
%cond al contorno
  VertexValue = [0 0];
  Vertex = [0 0];
  D = [0 0];
  InputVertexvalue = [0 0];
  
  %idxV = 1:length(geom.input.BC.InputVertexValues); %indice dei vertici

  % Se il lato e` di bordo
  if( any(e==-1) )
    % Lato di Dirichlet
	%in ogni nodo contiente il marker nodelist con la convezione che ha 0
    %se ha associato gdl altrimenti numero dispari se e di bordo che può
    %essere ereditato dal marker del lato che ho settato nella
    %triangolazione
   if( geom.pivot.nodelist(n(1))~=0 && geom.pivot.nodelist(n(2))~=0 ) %considero marker nodi estremi lato, se sono entrambi zero il lato è su neumann e tratto dopo
       %se uno dei marker e diverso da zero allora un lato è di neumann con
       %una cond di Dirichlet oppure il lato intero ha cond di Dirichlet e
       %devo imporre cond anche sul nodo medio
%-----------------------
   if( geom.pivot.nodelist(n(1)) ~= geom.pivot.nodelist(n(2)) )%se hanno lo stesso marker
       %devo cpaire che marker mettere nel punto medio
%
	if( any(geom.input.BC.InputVertexValues==geom.pivot.nodelist(n(1))) )% se il mio lato comincia in uno dei vertici del dominio
        %e ha 2 marker diversi e devo capire se è lato di neumann (un
        %vertice ha 0 e uno 5 ad esempio) 
        %se qualche marker degli estremi è uno dei nodi inziali
          % Vertice(1) con marker speciale
	  VertexValue(1) = 1; %se è così setto vertex=1 e so che mi trovo nei vertici del domino iniziale E MI SALVO 1 MARKE SPECIALE
          % e` il vettore con un 1 in corrispondenza del vertice
          % del lato con marker speciale
	  
	  Vertex(1) = geom.input.BC.Boundary.Values*...  %MARKER LATO CHE INIZIA NEL VERTICE 1 E VERIFICO CHE ABBAI IL MARKE DEL LATO A CUI è ASSOCIATO
              (geom.input.BC.InputVertexValues == geom.pivot.nodelist(n(1)))';%IL lato viene contraddistinto dal numero prmo vertice
          % valore del marker del lato del poligono iniziale che segue n1
          %devo fare attenzione all'uktimo lato
	  %InputVertexValue IN CUI SALVO IL MARKER DEL NODO
          InputVertexValue = geom.input.BC.InputVertexValues*...
	      (geom.input.BC.InputVertexValues == geom.pivot.nodelist(n(1)))';
          % MERKER nodo del poligono iniziale ce corrisponde al nodo n(1) del mio lato

	  % valore del marker del vertice del poligono iniziale n2
	  D(1) = geom.pivot.nodelist(n(2));%MARKER ASSOCIATO AL NODO 2 CHE è POTENZIALMENTE MARKER LATO
	end
%	
	if( any(geom.input.BC.InputVertexValues==geom.pivot.nodelist(n(2))) ) %NODO 2 è MARKR SPECIAL?
	  VertexValue(2) = 1; 
	  Vertex(2) = geom.input.BC.Boundary.Values*...
	  (geom.input.BC.InputVertexValues == geom.pivot.nodelist(n(2)))';
	  InputVertexValue = geom.input.BC.InputVertexValues*...
	      (geom.input.BC.InputVertexValues == geom.pivot.nodelist(n(2)))';
	  D(2) = geom.pivot.nodelist(n(1)); %SALVO MARKER NODO 1
   end
%VERIFICO I CASI DIVERSI CHE POSSONO VERIFICARSI
	if( sum(VertexValue) ~= 2 ) %IN TAL CASO NON HO 2 MARKER SPECIALI NEI NODI DEL LATO
	  Di = VertexValue*D'; %PRENDO IL MARKE DELL'ALTRO ESTREMO (NON QUELLO SPECIALE)
	  % nodo con condizione di Dirichlet
          geom.pivot.nodelist(nnode)= Di; %assegno al nuovo nodo il marke rdel lato 
          geom.pivot.Di(end+1,:) = [nnode, Di]; %INCREMENTO STRUTTURA INDICE NODO MARKE NODO
          geom.pivot.pivot(nnode) = min(geom.pivot.pivot)-1; %METTO UN NUMERO NEGATIVO MINORE DI TUTTI GLI ALTRI
	else
          % diamo al nuovo nodo il marker del lato
          % il lato che stiamo analizzando e` un lato del poligono
          % iniziale:
          
	  % l'indice del lato e` quello del nodo di inizio di quel lato
      %IN TAL CASO ENTRmambi gli extremi hano marker special
	  if( max(InputVertexvalue)-min(InputVertexvalue)>1 ) % siamo sul lato di chiusura IN INPUTVERTEXVALUE HO IL NUMERO DEL VERTICE ASSOCIATO AL VERTEX CON MARKER SPECIALE 1
	    %CONSIDERO IL MARKER DEL NODO PIù ALTO PER CONTRADD IL LATO
          Di = geom.input.BC.Boundary.Values(max(Vertex));%mi troov in caso in cui il triangolatore non ha diviso il lato
	  else % siamo sui lati 1->2->3->4->
	    Di = geom.input.BC.Boundary.Values(min(Vertex)); 
	  end
          % check della condizione di Neumann aperta
      if( rem(Di,2)== 0 ) % nodo con grado di liberta`, lato di
                              % Dirichlet aperto
	    geom.pivot.nodelist(nnode)= 0;
	    geom.pivot.pivot(nnode) = max(geom.pivot.pivot)+1;
      else
            geom.pivot.nodelist(nnode)= Di;
            geom.pivot.Di(end+1,:) = [nnode, Di];
            geom.pivot.pivot(nnode) = min(geom.pivot.pivot)-1;
      end
    end % if( sum(VertexValue) ~= 2 )
%----------------------------------
	else % if( geom.pivot.nodelist(n(1)) ~= geom.pivot.nodelist(n(2)) ) SE HANNO STRSSO MARKER ALLORA MI
        % TROVO SU UN LATO CHE STA SUL BORDO GENERATO DAL T RIANGOLATORE E POSSO USARE MARKER DI UNO DEI DUE ESTREMI
	     Di = geom.pivot.nodelist(n(1));
             geom.pivot.nodelist(nnode)= Di;
             geom.pivot.Di(end+1,:) = [nnode, Di];
             geom.pivot.pivot(nnode) = min(geom.pivot.pivot)-1;
             %HO AGGIORNATO LA STRUTTURE
    end % if( geom.pivot.nodelist(n(1)) ~= geom.pivot.nodelist(n(2)) )
%----------------------------------
%qui vado ad aggiornare le strutture nel caso di condizioni non omogenee e
%neumann
   else %SE ENTRO QUI SIGNIFICA CHE HO UN LATO CON VERTICE DIRICHLET E UNO NEUMANN
    % Lato di Neumann
      geom.pivot.nodelist(nnode) = 0;
      geom.pivot.pivot(nnode) = max(geom.pivot.pivot)+1;
    end % if( geom.pivot.nodelist(n(1))~=0 & geom.pivot.nodelist(n(2))~=0 )
      
    else % if( any(e==-1) ) %IN TAL CASO HO LATO INTERNO
      geom.pivot.nodelist(nnode) = 0;
      geom.pivot.pivot(nnode) = max(geom.pivot.pivot)+1;
  end %if( any(e==-1) )


end % for l=1:geom.nelements.nBorders