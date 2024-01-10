alpha = 10000;
beta = 4;
t = 1;
f = @(t, alpha, beta)  alpha/(1+exp(-beta*t));
df = @(t, alpha, beta) (alpha*beta*exp(-beta*t))/(exp(-beta*t) + 1)^2;
f(t,alpha,beta)
df(0.5, alpha, beta)