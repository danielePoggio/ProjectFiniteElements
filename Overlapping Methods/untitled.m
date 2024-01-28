% Definizione delle dimensioni
n = ...; % Dimensione della matrice A
m = ...; % Dimensione del blocco

% Definizione di A_{\Gamma \Gamma}, A_{\Gamma,1}, A_{1,1}, A_{1,\Gamma}, A_{\Gamma,2}, A_{2,2}, A_{2,\Gamma}
A_GammaGamma = ...; % Inserisci la tua implementazione
A_Gamma1 = ...;     % Inserisci la tua implementazione
A_11 = ...;         % Inserisci la tua implementazione
A_1Gamma = ...;     % Inserisci la tua implementazione
A_Gamma2 = ...;     % Inserisci la tua implementazione
A_22 = ...;         % Inserisci la tua implementazione
A_2Gamma = ...;     % Inserisci la tua implementazione

% Creazione della funzione di prodotto matrice-vettore
A_product = @(x) A_GammaGamma*x - A_Gamma1*(A_11\x) - A_Gamma2*(A_22\x);

% Creazione della funzione di precondizionamento basata su Cholesky incompleta
preconditioner = ichol(A_11, struct('michol','on')) + ichol(A_22, struct('michol','on'));

% Definizione del vettore destro b
b = ...; % Inserisci il tuo vettore destro

% Applicazione del PCG con precondizionatore
[x, flag, relres, iter, resvec] = pcg(A_product, b, 1e-6, 100, preconditioner);

% Visualizzazione dei risultati
disp(['Flag di convergenza: ', num2str(flag)]);
disp(['Numero di iterazioni: ', num2str(iter)]);
disp(['Residuo relativo finale: ', num2str(relres)]);

% Plot del residuo nel corso delle iterazioni
figure;
semilogy(resvec, '-o');
title('Convergenza del PCG');
xlabel('Iterazione');
ylabel('Residuo relativo');
