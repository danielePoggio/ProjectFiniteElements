function y = funzione_specifica(x)
    % Funzione monotona crescente nell'intervallo [0, 1]
    % con derivata prima massima in x = 0.5

   % Definizione della funzione simile ad una sigmoide
    a = 10;  % Parametro di regolazione della "steepness"
    y = 1 ./ (1 + exp(-a * (x - 0.5)));
end


