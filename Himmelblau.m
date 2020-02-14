function cost = Himmelblau(z)
    x = squeeze(z(1,:,:));
    y = squeeze(z(2,:,:));
   cost = (x.^2 +y-11).^2 + (x+y.^2-7).^2;
end