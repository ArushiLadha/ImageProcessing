function approx = escapeVelocity(z0, c, N)
    next = 0;
    pre = z0;
    approx = z0;
    n = 1;
    while abs(pre)<2 && n<N
        next = pre^2 + c;
        %next = exp(z^3) - c;
        n = n+1;
        pre = next;
        approx = log(abs(next))/(2^n);
    end;
end

%exp(z3) - 0.621