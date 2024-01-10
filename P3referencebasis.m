base = @(x,y) [x^3 y^3 (x^2)*y x*(y^2) x^2 y^2 x*y x y 1]';
loc = @(x,y) [x^3 y^3 (x^2)*y x*(y^2) x^2 y^2 x*y x y 1;
    3*x^2 0 2*x*y y^2 2*x 0 y 1 0 0;
    0 3*y^2 x^2 2*x*y 0 2*y x 0 1 0];
matrix = [loc(0,0);
    loc(1,0);
    loc(0,1);
    base(1/3,1/3)'];

weights = zeros(10,10);
b = eye(10,10);
for i=1:10
    w = matrix\b(:,i);
    weights(i,:) = w';
end

phiP3 = @(x,y) weights*base(x,y);


run("")