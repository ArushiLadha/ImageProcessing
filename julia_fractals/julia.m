function M = julia(zmax,c,N)
    x = linspace(-zmax, zmax,500);
    y = linspace(-1i*zmax, 1i*zmax,500);
    [X, Y] = meshgrid(x,y);
    Z = Y+X;
    M = zeros(500);
    for a = 1:500
        for b = 1:500
            M(a,b) = escapeVelocity(Z(a,b), c, N);            
        end
    end
    imagesc(atan(0.1*M));
    axis xy;
end
    