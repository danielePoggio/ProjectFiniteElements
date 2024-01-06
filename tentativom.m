% Genera una griglia di punti nel quadrato [0,1]x[0,1]
x = linspace(0, 1, 100);
y = linspace(0, 1, 100);
[X, Y] = meshgrid(x, y);

% Calcola i valori della funzione su tutta la griglia
Z = funzione_specifica(X, Y);

% Visualizza la funzione tridimensionale
figure;
surf(X, Y, Z, 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('Valore');
title('Funzione con Gradiente Pendente nel Mezzo');

% Mostra il grafico