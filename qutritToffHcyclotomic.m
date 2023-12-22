%% Implementation of qutrit controlled Clifford+T gates
% For the paper Exact synthesis of multiqutrit Clifford-cyclotomic circuits
% An earlier implementation, with fewer gates and where some of the qutrit
% Clifford gates differ in definition by a global phase, is at https://git
% hub.com/lia-approves/qudit-circuits/tree/main/qutrit_control_Clifford_T

clear; % clear all variables so only the minimal gate set is in scope
global d I II X H CCX; % declare the minimal gate set

d = 3;
I = eye(d);
II = kron(I,I);
X = init_X();
H = init_H();
CCX = init_CCX();

ZCXupsidedown = init_ZCX(true);
CXupsidedown = multAll(arrayfun(@(j) kron(I,X^j) * ZCXupsidedown^j * kron(I,X^(d-1)^j), 1:d-1, 'UniformOutput', false));
swap = kron(I,H^2)*CX * kron(H^2,I)*CXupsidedown * kron(I,H^2)*CX;

P3 = swap*ZCX*swap * ZCX * swap*ZCX^(d-1)*swap * ZCX^(d-1);
Z_OCX01 = P3*swap*CX*swap*P3*swap*CX^(d-1)*swap;
X01 = kron(k0',I)*Z_OCX01*kron(k0,I);
ZCX01 = kron(X,X01) * (Z_OCX01* kron(X^2,I))^((d-1)/2);
CX01 = kron(X,I) * (ZCX01*kron(X^2,I))^((d-1)/2);
ZZCX01 = kron(I,CX01) * kron(ZCX,I) * kron(I,CX01) * kron(ZCX^(d-1),I) * kron(swap,I)*kron(I,ZCX01)*kron(swap,I);
ZZCX = (ZZCX01 * kron(II,X))^(d-1) * kron(II,X);

Q0 = kron(I,k0'*X^(d-1)*H) * ZCX * kron(I,H^3*X*k0);

%if d = 3
    ZZCw = (kron(II,H^3)*ZZCX*kron(II,H) * kron(II,X^(d-1) * X01 * X))^2;
    %ZZCwdag = (kron(II,H)*ZZCX*kron(II,H^3) * kron(II,X^(d-1) * X01 * X))^2;
    ZCS = removeBorrowedAncilla(kron(kron(I,X^(d-1)),I) * ZZCw * kron(kron(I,X),I),1);
    TCiHdag = kron(X^(d-1),H^2) * ZCS * kron(I,H^2) * ZCS * kron(I,H^3) * ZCS * kron(I,H^2) * ZCS * kron(I,H) * ZCS * kron(I,H^2) * ZCS * kron(X,I);
    %TCminusiH = kron(X^(d-1),H^2) * ZCS^2 * kron(I,H^2) * ZCS^2 * kron(I,H^3) * ZCS^2 * kron(I,H^2) * ZCS^2 * kron(I,H) * ZCS^2 * kron(I,H^2) * ZCS^2 * kron(X,I);
    RtensI = TCiHdag^2 * kron(X^(d-1),X)*ZCX01*kron(X,X^(d-1));
    R = removeBorrowedAncilla(RtensI,1);
%end

1;

function X = init_X()
    global I;
    X = I(:,[2:end,1]);
end

function H = init_H()
    global d;
    w = exp(1i*2*pi/d);
    Hglobalphase = -w^2/sqrt(-d); % This global phase is irrelevant to constructing R, but affects the global phase of other constructed gates such as Hdag
    H = Hglobalphase * ( w .^ mod( arrayfun(@(x,y) x*y, repmat((0:d-1).',1,d), repmat(0:d-1,d,1)), d));
end

function CCX = init_CCX()
    global d X;
    CX = cellfun(@(k) X^k,num2cell([0:d-1]),'UniformOutput',false);
    CX = blkdiag(CX{:});
    CCX = cellfun(@(k) CX^k,num2cell([0:d-1]),'UniformOutput',false);
    CCX = blkdiag(CCX{:});
end

function Hdag = Hdag()
    global d H;
    Hdag = H^(4*d-1); % This is for exact global phase. Up to a global phase, H^3 suffices
end

function H2m = H2m() % H^2 but with a global phase of -1. If the global phase were 1, sends x mod d to -x mod d.
    global d H;
    H2m = H^(2*d);
end

function CX = CX()
    global d II X H CCX;
    CX = (kron(X * H2m() * X^(d-1), II) * CCX^((d+1)/2))^2;
    CX = removeBorrowedAncilla(CX,3);
end

function CXdag = CXdag()
    global d;
    CXdag = CX()^(d-1);
end

function ZCX = ZCX()
    global d I X H CCX;
    ZCX = kron(kron(I,Hdag()^2),I) * CCX * kron(kron(I,H^2),I) * kron(CXdag(),I) * CCX * kron(CX(),I);
    ZCX = removeBorrowedAncilla(ZCX,2);
    ZCX = kron(I,X) * ZCX^(d-1);
end

function V = removeBorrowedAncilla(U,b) % U is the unitary, and b is the index of the borrowed ancilla
    global d I;
    idxs = arrayfun(@(k) k+1:k+d^(b-1), 0:d^b:size(U,2)-1, 'UniformOutput', false);
    V = U([idxs{:}],[idxs{:}]);
    Uv = arrayfun(@(x,y) U([idxs{:}]+x*d^(b-1),[idxs{:}]+y*d^(b-1)), 0:d-1, 0:d-1, 'UniformOutput', false);
    if any(arrayfun(@(k) any(any(V-[Uv{k}] > 1e-5)), 1:d))
        strV = sprintf([repmat('%.3f ', 1, size(V, 2)) '\n'], V);
        strU = sprintf([repmat('%.3f ', 1, size(U, 2)) '\n'], U);
        ME = MException('component:noBorrowedAncilla',[strV 'not a borrowed ancilla of ' strU], strV, strU);
        throw(ME)
    end
end

function prod = multAll(list)
    prod = list{1};
    for i = 2:numel(list)
        prod = list{i} * prod;
    end
end
