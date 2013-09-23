function M = juliaAnimation(zmax,c,N)
    x = linspace(-zmax, zmax,500);
    y = linspace(-1i*zmax, 1i*zmax,500);
    [X, Y] = meshgrid(x,y);
    Z = Y+X;
    M = N*ones(500);
    clear x y X Y;
    pre = Z;
    counter = 0;
    C = c*ones(500);
    close all;
%     colormap(hot);
    h = imagesc(0*M);
    pause(.1);
    for n = 1:N
        for a = 1:500
            for b = 1:500
                if M(a,b)==N && abs(pre(a,b))>3
                    M(a,b) = n + 1 - log(log(abs(pre(a,b))))/log(2);
                    counter = counter+1;
                end
                if counter == 250000
                    a = 500; n =N;
                    break;
                end
            end
        end
        set(h,'CData',atan(0.1*M));
        drawnow;
        pre = pre.^2 + C;       
    end
    
    
    