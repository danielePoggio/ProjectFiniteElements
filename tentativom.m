% f = @(x,K,grade) K*x.^grade./(K+x.^grade);
% Genera una griglia di punti nel quadrato [0,1]x[0,1]
x = linspace(0, 1, 100);
y = linspace(0, 1, 100);
[X, Y] = meshgrid(x, y);
f = @(x,y) funzione_specifica(x+y);
% Calcola i valori della funzione su tutta la griglia
Z = f(X,Y);
figure;
surf(X, Y, Z, 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('Valore');
title('Funzione con Gradiente Lungo la Diagonale y=x');