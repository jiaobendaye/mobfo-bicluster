x = linspace(-3, 3);
y = linspace(-3, 3);
[X, Y] = meshgrid(x, y);

Z = sin(X) + cos(Y);
% Z = GoldsteinPrice([X;Y]);
figure;
contour(X,Y,Z)
