clear all
close all
clc

%% Eseguo Triangolazione sul Dominio
area = 0.02;
geom = Triangolator(area);
close all

%% Problema differenziale
u0 = @(x,y) 16*x*(1-x)*y*(1-y);
gradu = @(x,y) 16*[y*(1-y)*(1 -2*x), x*(1-x)*(1 - 2*y)]';
rho = @(x,y) 1.0;
mu = @(x,y) 1.0;
beta = @(x,y) [0.0, 0.0];
sigma = @(x,y) 0.0;
f = @(t,x,y) 32*(x*(1-x) + y*(1-y)) + 16*(1-2*x)*(y*(1-y));
gNe = @(t,x,y) -16*x*(1-x);
gDi = @(t,x,y) 0;
T = 5;
dt = 0.1;
Nt = fix(T/dt); % parte intera 


%% Soluzione problema discretizzato
uh = euleroImplicit(geom, rho, mu, beta, sigma, f, gDi, gNe, u0, dt,Nt);
