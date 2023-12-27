%% Implementation of qutrit Toffoli+Hadamard gate set level operators
% In paper: Exact synthesis of multiqutrit Clifford-cyclotomic circuits
% For each level operator (-1)_[a], \omega_[a], X_[a,b], H_[a,b,c]: Add
% controls or conjugate by classical reversible gates constructions here
% for exact implementation on any number of qutrits + 2 borrowed ancilla

clear; % clear all variables so only the minimal gate set is in scope
global d I II R1 H CCX; % declare the minimal gate set

d = 3;
I = eye(d);
II = kron(I,I);
R1 = init_R1();
H = init_H();
CCX = init_CCX();

% Step through each line below to check its truth table (for gates that
% are a diagonal gate times a d-ary classical reversible gate) or matrix
truthU(R1);
disp(H);
truthU(CCX);
disp(Hdag());
truthU(X());
truthU(Xdag());
truthU(Hsqm());
truthU(wgp());
truthU(CX());
truthU(CXdag());
truthU(Swapm());

% The constructions from here on are specific to d = 3. For d > 3 it's
% easier to use the gate set {X, H, ZCX} instead of {X, H, CCX}.
if d == 3
    truthU(ZCX());
    truthU(ZCXdag());
    truthU(S());
    truthU(mgp());
    truthU(Swap());
    truthU(X0_1());
    truthU(X00_01());
    truthU(X000_001());
    truthU(w22());
    disp(Hdag2mwsq());
    truthU(M2());
    disp(H2());
    disp(Hdag2());
end

%% Except for init functions, all functions here are constructions using
% only the constructions in functions above it or in the same folder as it.

function R1 = init_R1()
    global d;
    R1 = diag(arrayfun(@(k) exp(1i*2*pi/d)^k, 0:d-1));
end

function H = init_H()
    global d;
    w = exp(1i*2*pi/d);
    Hglobalphase = -w^2/sqrt(-d); % This global phase is irrelevant to
    % building (-1)_[2], but affects global phase of other constructions
    H = Hglobalphase * ( w .^ mod( arrayfun(@(x,y) x*y, ...
        repmat((0:d-1).',1,d), repmat(0:d-1,d,1)), d));
end

function CCX = init_CCX()
    global d I;
    X = I(:,[2:end,1]);
    CX = cellfun(@(k) X^k,num2cell([0:d-1]),'UniformOutput',false);
    CX = blkdiag(CX{:});
    CCX = cellfun(@(k) CX^k,num2cell([0:d-1]),'UniformOutput',false);
    CCX = blkdiag(CCX{:});
end

% Hdag with exact global phase. Is H^3 up to a global phase.
function Hdag = Hdag()
    global d H;
    Hdag = H^(4*d-1);
end

function X = X()
    global R1 H;
    X = Hdag() * R1 * H;
end

function Xdag = Xdag()
    global d;
    Xdag = X()^(d-1);
end

% Sends |x> to -|-x mod d>. Is H^2 up to a global phase.
function Hsqm = Hsqm()
    global d H;
    Hsqm = H^(2*d);
end

function wgp = wgp()
    global d H;
    [~,C,~] = gcd(8,d);
    wgp = H^(4*mod(C,d));
end

function CX = CX()
    global d II H CCX;
    CX = (kron(X() * Hsqm() * Xdag(), II) * CCX^((d+1)/2))^2;
    CX = removeBorrowedAncilla(CX,1);
end

function CXdag = CXdag()
    global d;
    CXdag = CX()^(d-1);
end

% Swap gate, with a global phase of -1.
function Swapm = Swapm()
    global I;
    udCX = reIndex(CX(), [2 1]); % upside-down CX gate
    Swapm = kron(I,Hsqm())*CX() *kron(Hsqm(),I)*udCX* kron(I,Hsqm())*CX();
end

%% All below here are for d = 3 and do not necessarily hold for d > 3.

% |0>-controlled X gate: Does to target X if control is |0>, else does I
function ZCX = ZCX()
    global d I H CCX;
    nZCX = kron(kron(I,Hdag()^2),I) * CCX * kron(kron(I,H^2),I) * ...
        kron(CXdag(),I) * CCX * kron(CX(),I); % |not 0>-controlled X
    nZCX = removeBorrowedAncilla(nZCX,2);
    ZCX = kron(I,X()) * nZCX^(d-1);
end

% |0>-controlled Xdag: Does to target Xdag if control is |0>, else I
function ZCXdag = ZCXdag()
    global d;
    ZCXdag = ZCX()^(d-1);
end

% The 1-qutrit S gate. Is the \omega_[2] operator applying \omega to |2>
function S = S()
    global H;
    S = (kron(X()^2,H)*ZCX()*kron(X(),H^3*X()^2*H^2*X()))^2;
    S = removeBorrowedAncilla(S,2);
end

% -1 global phase 1-qutrit gate
function mgp = mgp()
    global d H;
    mgp = (S()*H)^(d^2);
end

% Swap 2-qutrit gate
function Swap = Swap()
    global I;
    Swap = Swapm() * kron(mgp(),I);
end

% X_[0,1] level operator which is a 1-qutrit gate, for d = 3.
function X0_1 = X0_1()
    X0_1 = mgp()*X()^2*Hsqm()*X();
end

% X_[00,01] level operator, i.e. |0>-controlled X_[0,1] gate, for d = 3
% which is unitary.  For d > 3, known constructions have computational
% basis initialization and deterministic measurement due to X_[0,1].
function X00_01 = X00_01()
    global d I H;
    udZCX = reIndex(ZCX(),[2 1]); % upside-down ZCX gate
    udZCXdag = reIndex(ZCXdag(),[2 1]);
    perm00_10_01 = udZCX * ZCX() * udZCXdag * ZCXdag();
    udCX = reIndex(CX(),[2 1]);
    udCXdag = reIndex(CXdag(),[2 1]);
    X00_01 = kron(X(),X0_1()) * (perm00_10_01 * udCX * ...
        perm00_10_01 * udCXdag * kron(X()^2,I))^((d-1)/2);
end

function X000_001 = X000_001()
    global d I;
    CX0_1 = (kron(Xdag()^2,I) * X00_01)^((d-1)/2) * kron(Xdag(),I);
    X000_001 = kron(I,CX0_1) * kron(ZCX(),I) * kron(I,CX0_1) * ...
        kron(ZCXdag(),I) * reIndex(X00_01(),[1 3 2]);
end

% The \omega_[22] level operator which applies phase \omega to |22>.
% Is the |2>-controlled S gate.
function w22 = w22()
    global d II H;
    ZZCX = kron(II,X())*(kron(II,X())*X000_001())^(d-1); % |00>-ctrled X
    w22 = ( kron(kron(X()^2,X()^2),H) * ZZCX * ...
        kron(kron(X(),X()),H^3*X()^2*H^2*X()) )^2;
    w22 = removeBorrowedAncilla(w22,3);
end

% The Hdag_[20,21,22] level operator up to controlled global phase -w^2
function Hdag2mwsq = Hdag2mwsq()
    global I H;
    Zwsqwsq = kron(I,H^2) * w22() * kron(I,H^2) * w22();
    Hdag2mwsq = Zwsqwsq * kron(I,Hdag())*Zwsqwsq*kron(I,H) * Zwsqwsq;
end

% The (-1)_[2] level operator
function M2 = M2()
    M2 = kron(Xdag(),X())*X00_01()*kron(X(),Xdag()) * Hdag2mwsq()^2;
    M2 = removeBorrowedAncilla(M2,2);
end

% The H_[20,21,22] 2-qutrit level operator
function H2 = H2()
    global I H;
    Zww = w22()^2 * kron(I,Hdag()^2) * w22()^2 * kron(I,Hdag()^2);
    H2mw = Zww * kron(I,Hdag()) * Zww * kron(I,H) * Zww;
    H2 = H2mw * kron(M2() * S()^2, I);
end

% The Hdag_[20,21,22] 2-qutrit level operator
function Hdag2 = Hdag2()
    global I;
    Hdag2 = Hdag2mwsq() * kron(M2() * S(), I);
end