function valore = funzione_specifica(x, y)
    % Calcola il gradiente lungo la diagonale y=x
    gradiente_diagonale = 50;

    % Calcola la distanza dalla diagonale
    distanza_diagonale = abs(y - x);

    % Calcola la funzione
    valore = exp(-gradiente_diagonale * (distanza_diagonale).^2);
end

