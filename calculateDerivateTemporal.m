% Step 1: Declare symbolic variable
syms t x y;

% Step 2: Create a symbolic expression

% Step 3: Differentiate the expression
expr_diff_t = diff(u, t);
expr_diff_x = diff(u, x);
expr_diff_y = diff(u, y);

% Step 4: Convert the symbolic expression to a string
expr_diff_t_str = char(expr_diff_t);
expr_diff_x_str = char(expr_diff_x);
expr_diff_y_str = char(expr_diff_y);

% Step 5: Create an anonymous function using eval and str2fun
dut = str2func(['@(t,x,y)'  expr_diff_t_str ]);
gradu = str2func(['@(t,x, y) ' '[(' expr_diff_x_str '), (' expr_diff_y_str ')]']);

% Ripetiamo per creare l'Hessiana
expr_diff_x = diff(gradu, x);
expr_diff_y = diff(gradu, y);
expr_diff_x_str = char(expr_diff_x);
expr_diff_y_str = char(expr_diff_y);
Hu = str2func([ '@(t,x,y)'  '[' expr_diff_x_str '  ;  ' expr_diff_y_str ']']);

clear dux duy d2uy d2ux expr_diff_y_str expr_diff_x_str expr_diff_t_str expr_diff_y expr_diff_x expr_diff_t graducell t x y
