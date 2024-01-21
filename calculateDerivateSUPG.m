clear all
clc
syms x y;

N1 = @(x,y) x;
N2 = @(x,y) y;
N3 = @(x,y) 1 - x - y;

phi1 = @(x,y) 2*N1(x,y)*(N1(x,y) - 0.5);
phi2 = @(x,y) 2*N2(x,y)*(N2(x,y) - 0.5);
phi3 = @(x,y) 2*N3(x,y)*(N3(x,y) - 0.5);
phi4 = @(x,y) 4*N3(x,y)*N1(x,y);
phi5 = @(x,y) 4*N1(x,y)*N2(x,y);
phi6 = @(x,y) 4*N2(x,y)*N3(x,y);
phi = @(x,y) [phi1(x,y), phi2(x,y), phi3(x,y), phi4(x,y), phi5(x,y), phi6(x,y)]';
dxphi_syms = diff(phi,x);
expr_dxphi = char(dxphi_syms);
dxphi = str2func(['@(x,y)' expr_dxphi]);

dyphi_syms = diff(phi,y);
expr_dyphi = char(dyphi_syms);
dyphi = str2func(['@(x,y)' expr_dyphi]);

Jphi = @(x,y) [dxphi(x,y) dyphi(x,y)];

d2xphi_syms = diff(Jphi, x);
expr_d2xphi = char(d2xphi_syms);
d2xphi = str2func(['@(x,y)' expr_d2xphi]);

d2yphi_syms = diff(Jphi,y);
expr_d2yphi = char(d2yphi_syms);
d2yphi = str2func(['@(x,y)' expr_d2yphi]);