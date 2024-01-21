function [gradu, d2u] = calculateDerivate(u)
% Step 1: Declare symbolic variable
syms x y;

% Step 2: Create a symbolic expression

% Step 3: Differentiate the expression
expr_diff_x = diff(u, x);
expr_diff_x2 = diff(expr_diff_x, x);
expr_diff_y = diff(u, y);
expr_diff_y2 = diff(expr_diff_y, y);

% Step 4: Convert the symbolic expression to a string
expr_diff_x_str = char(expr_diff_x);
expr_diff_y_str = char(expr_diff_y);

% Step 5: Create an anonymous function using eval and str2fun
gradu = str2func(['@(x, y) ' '[' expr_diff_x_str '; ' expr_diff_y_str ']']);

% Ripetiamo per creare l'Hessiana
expr_diff_x_str = char(expr_diff_x2);
expr_diff_y_str = char(expr_diff_y2);
d2u = str2func([ '@(x,y)'   expr_diff_x_str '  +  ' expr_diff_y_str ]);

end