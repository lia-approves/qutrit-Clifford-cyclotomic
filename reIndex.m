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