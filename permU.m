% Returns permutation gate with old qudit indices mapped to the indices in
% maplist.  Only use this function after having the 2-qudit swap gate,
% with the exception of for re-indexing a gate.
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