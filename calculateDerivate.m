% Step 1: Declare symbolic variable
syms x y t;

% Step 2: Create a symbolic expression

% Step 3: Differentiate the expression
expr_diff_t = diff(u, t);
expr_diff_x = diff(u, x);
expr_diff_y = diff(u, y);

% Step 4: Convert the symbolic expression to a string
expr_diff_t_str = char(expr_diff_t);
expr_diff_x_str = char(expr_diff_x);
expr_diff_y_str = char(expr_diff_y);

% Step 5: Create an anonymous function using eval and str2func
dut = str2func(['@(t,x,y) ' expr_diff_t_str]);
dux = str2func(['@(t,x,y) ' expr_diff_x_str]);
duy = str2func(['@(t,x,y) ' expr_diff_y_str]);
gradu = @(t,x,y) [duy(t,x,y), dux(t,x,y)]';
graducell = @(t,x,y) {duy, dux}';

% Ripetiamo per creare l'Hessiana
expr_diff_x = diff(graducell, x);
expr_diff_y = diff(graducell, y);
expr_diff_x_str = char(expr_diff_x);
expr_diff_y_str = char(expr_diff_y);
d2ux = str2func(['@(t,x,y) ' expr_diff_x_str]);
d2uy = str2func(['@(t,x,y) ' expr_diff_y_str]);
Hu = @(t,x,y) [d2ux(t,x,y), d2uy(t,x,y)];

clear dux duy d2uy d2ux expr_diff_y_str expr_diff_x_str expr_diff_t_str expr_diff_y expr_diff_x expr_diff_t graducell t x y

