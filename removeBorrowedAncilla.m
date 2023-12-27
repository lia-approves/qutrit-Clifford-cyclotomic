% Returns unitary U except without the borrowed ancilla with index b.
% If the qudit with index b is not a borrowed ancilla, throws error.
function V = removeBorrowedAncilla(U,b)
    global d I;
    b = logd(size(U,1))-b+1;
    idxs = arrayfun(@(k) k+1:k+d^(b-1), 0:d^b:size(U,2)-1, ...
        'UniformOutput', false);
    V = U([idxs{:}],[idxs{:}]);
    Uv = arrayfun(@(x,y) U([idxs{:}]+x*d^(b-1),[idxs{:}]+y*d^(b-1)), ...
        0:d-1, 0:d-1, 'UniformOutput', false);
    if any(arrayfun(@(k) any(any(V-[Uv{k}] > 1e-5)), 1:d))
        strV = sprintf([repmat('%.3f ', 1, size(V, 2)) '\n'], V);
        strU = sprintf([repmat('%.3f ', 1, size(U, 2)) '\n'], U);
        ME = MException('component:noBorrowedAncilla',[strV ['not a' ...
            'borrowed ancilla of '] strU], strV, strU);
        throw(ME)
    end
end