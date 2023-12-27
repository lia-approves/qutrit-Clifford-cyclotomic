% Return 2-qudit swap gate on n qudits. This function is only used after
% having the 2-qudit swap gate, with the exception of re-indexing a gate.
function swapn = swapU(i,j,n)
    global d I;
    if j < i
        swapn = swapU(j,i,n);
    elseif i == j
        swapn = eye(d^n);
    else
        swapn = zeros(d^n);
        for a = 1:d
            for b = 1:d
                swapn = swapn + tensall({tenpow(I,i-1),I(:,a),tenpow(I, ...
                    j-i-1),I(:,b),tenpow(I,n-j)}) * tensall({tenpow(I, ...
                    i-1),I(:,b),tenpow(I,j-i-1),I(:,a),tenpow(I,n-j)})';
            end
        end
    end
end