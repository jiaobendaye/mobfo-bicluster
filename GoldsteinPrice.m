function cost = GoldsteinPrice(z)
    x = squeeze(z(1,:,:));
    y = squeeze(z(2,:,:));
    cost = (1+((x+y+1).^2).*(19-14*x+3*x.^2-14*y+6*x.*y + 3*y.^2)).*...
        (30+((2*x-3*y).^2).*(18-32*x+12*x.^2+48*y-36*x.*y+27*y.^2));
end