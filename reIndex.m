% Returns U conjugated by permutation gates, with old qudit indices mapped
% to the indices in maplist.
function newU = reIndex(U,maplist)
    global d;
    n = length(maplist);
    if size(U,1) < d^n
        U = kron(U, eye(d^n / size(U,1)));
    end
    perm_gate = permU(maplist);
    newU = perm_gate * U * perm_gate;
end

% Returns permutation gate with old qudit indices mapped to the indices in
% maplist.  Only use this function after having the 2-qudit swap gate,
% with the exception of re-indexing a gate.
function perm = permU(maplist)
    global d;
    perm = eye(d^length(maplist));
    for a = 1:length(maplist)
        if a ~= maplist(a)
            perm = swapU(a,maplist(a),length(maplist))*perm;
            maplist([a maplist(a)]) = maplist([maplist(a) a]);
        end
    end
end

% Return 2-qudit swap gate on n qudits. Only use this function after
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