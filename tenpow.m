% tensor to the power.  returns base^{\otimes n}
function powd = tenpow(base,n)
    if n == 0
        powd = 1;
        return;
    end
    powd = 1;
    for i = 1:n
        powd = kron(base,powd);
    end
end