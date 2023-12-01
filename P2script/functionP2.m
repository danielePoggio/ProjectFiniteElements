%% Scriviamo le funzioni base di P2
Nv = 6;
N1 = @(x,y) x;
N2 = @(x,y) y;
N3 = @(x,y) 1 - x - y;
N = @(x,y) [N1(x,y), N2(x,y), N3(x,y)];

phi1 = @(x,y) 2*N1(x,y)*(N1(x,y) - 0.5);
phi2 = @(x,y) 2*N2(x,y)*(N2(x,y) - 0.5);
phi3 = @(x,y) 2*N3(x,y)*(N3(x,y) - 0.5);
phi4 = @(x,y) 4*N3(x,y)*N1(x,y);
phi5 = @(x,y) 4*N1(x,y)*N2(x,y);
phi6 = @(x,y) 4*N2(x,y)*N3(x,y);
phi = @(x,y) [phi1(x,y), phi2(x,y), phi3(x,y), phi4(x,y), phi5(x,y), phi6(x,y)];

Jphi = @(x,y) [4*x - 1, 0, 4*x + 4*y - 3, 4 - 4*y - 8*x, 4*y, -4*y; 
    0, 4*y - 1, 4*x + 4*y - 3, -4*x, 4*x, 4 - 8*y - 4*x]';

clear N1 N2 N3 phi1 phi2 phi3 phi4 phi5 phi6





% idxPhiFunction = zeros(dim,3); % vettore che mi serve per implementare le funzioni
% idx = 1;
% for i=1:2
%     for j=1:2
%         for k=1:2
%             if (i+j+k) == 2
%                 idxPhiFunction(idx,:) = [i, j, k];
%                 idx = idx + 1;
%             end
%         end
%     end
% end
% virtualIdx = idxPhiFunction; % creo una coppia su cui fare manipolazioni
