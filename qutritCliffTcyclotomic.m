%% Implementation of qutrit Clifford+T gate set level operators
% From paper: "Exact synthesis of multiqutrit Clifford-cyclotomic circuits"

% Any level operator (-1)_[a], \zeta_[a], X_[a,b], H_[a,b,c] has exact
% implementation here as a unitary over {H, CX, R2} + 2 borrowed ancilla,
% up to adding controls and/or conjugation by classical reversible gates
% ancilla-free in "Constructing All Qutrit Controlled Clifford+T gates in
% Clifford+T" doi:10.1007/978-3-031-09005-9_3 and "The Qudit ZH-Calculus:
% Generalised Toffoli+Hadamard and Universality" doi:10.4204/EPTCS.384.9

clear; % clear all variables so only the minimal gate set is in scope
global d I II H CX R2; % declare the minimal gate set

d = 3;
I = eye(d);
II = kron(I,I);
H = init_H();
CX = init_CX();
R2 = init_R2();

% Step through each line below to check its truth table (for gates that
% are a diagonal gate times a d-ary classical reversible gate) or matrix
disp(H); % The 1-qutrit H_[0,1,2] level operator
truthU(CX);
truthU(R2);
disp(Hdag());
truthU(CXdag());
truthU(R2dag());
truthU(R1());
truthU(X());
truthU(Xdag());
truthU(Hsqm());
truthU(wgp());
truthU(Swapm());
truthU(P1zdg(0));
truthU(ZCXutgp());
truthU(ZCP1dag0());

% The constructions from here on are specific to d = 3. For |0>-controlled
% X and controlled H in odd prime qudit Clifford+R2, see "Scaling W state
% circuits in the qudit Clifford hierarchy" doi:10.1145/3594671.3594687
if d == 3
    truthU(mgp()); % The (-1) global phase 1-qutrit operator
    truthU(Swap());
    truthU(X0_1()); % The X_[0,1] 1-qutrit level operator
    truthU(S());
    truthU(Sdag());
    truthU(S8());
    truthU(Sdag1());
    truthU(zgp()); % The \zeta global phase 1-qutrit operator
    truthU(zdaggp());
    truthU(ZCX);
    truthU(X00_01()); % The X_[00,01] 2-qutrit level operator
    truthU(X000_001()); % The X_[000,001] 3-qutrit level operator
    truthU(w22());
    truthU(wdag1_22());
    truthU(z2()); % The \zeta_[2] 1-qutrit level operator
    disp(H2m());
    disp(Hdag2m());
    truthU(m2()); % The (-1)_[2] 1-qutrit level operator
    disp(H2()); % The H_[20,21,22] 2-qutrit level operator
    disp(Hdag2());
end

%% Except for init functions, all functions here are constructions using
% only the constructions in functions above it or in the same folder as it.

function H = init_H()
    global d;
    w = exp(1i*2*pi/d);
    Hglobalphase = -w^2/sqrt(-d); % This global phase is irrelevant to
    % building (-1)_[2], but affects global phase of other constructions
    H = Hglobalphase * ( w .^ mod( arrayfun(@(x,y) x*y, ...
        repmat((0:d-1).',1,d), repmat(0:d-1,d,1)), d));
end

function CX = init_CX()
    global d I;
    X = I(:,[2:end,1]);
    CX = cellfun(@(k) X^k,num2cell([0:d-1]),'UniformOutput',false);
    CX = blkdiag(CX{:});
end

function R2 = init_R2()
    global d;
    R2 = diag(arrayfun(@(k) exp(1i*2*pi/d^2)^k, 0:d-1));
end

% Hdag with exact global phase. Is H^3 up to a global phase.
function Hdag = Hdag()
    global d H;
    Hdag = H^(4*d-1);
end

function CXdag = CXdag()
    global d CX;
    CXdag = CX()^(d-1);
end

function R2dag = R2dag()
    global d R2;
    R2dag = R2^(d^2-1);
end

function R1 = R1()
    global d R2;
    R1 = R2^d;
end

function X = X()
    global H;
    X = Hdag() * R1() * H;
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

% Swap gate, with a global phase of -1.
function Swapm = Swapm()
    global I CX;
    udCX = reIndex(CX, [2 1]); % upside-down CX gate
    Swapm = kron(I,Hsqm())*CX * kron(Hsqm(),I)*udCX * kron(I,Hsqm())*CX;
end

% 1-qudit P1(a) gate, with global phase \zetadag
function P1zdg = P1zdg(a)
    global R2;
    P1zdg = X()^a * R2dag() * X() * R2 * Xdag()^(a+1);
end

function ZCXutgp = ZCXutgp()
    global d I H CX R2;
    ZCXutgp = kron(P1zdg(0)^((d-1)/2),Hdag()) * (kron(I,R2)*CX)^d * kron(I,H);
end

% |0>-ctrled P1dag(0) up to controlled global phase and global phase.
function ZCP1dag0 = ZCP1dag0()
    global d I R2;
    ZCP1dag0 = ZCXutgp()^(d-1) * kron(I,R2) * ZCXutgp() * kron(I,R2dag());
    if d == 3 % then can correct controlled global phase which is Clifford
        ZCP1dag0 = ZCP1dag0 * kron(wgp(),I);
    end
end

%% All below here are for d = 3 and do not necessarily hold for d > 3.

function mgp = mgp()
    global d H;
    mgp = (P1zdg(2)*H)^(d^2);
end

% Swap 2-qutrit gate
function Swap = Swap()
    global I;
    Swap = Swapm() * kron(mgp(),I);
end

% 1-qutrit X_[a,b]
function X0_1 = X0_1(a,b)
    global H;
    if nargin < 2
        a = 0;
        b = 1;
    end
    if any(mod(a-b,3) == [1,2])
        X0_1 = mgp() * X()^(a+b)^2 * Hsqm() * X()^(a+b);
    else
        ME = MException('component:invalidAB', ...
            sprintf('Invalid a=%d, b=%d',a,b));
        throw(ME)
    end
end

function S = S()
    global I H CX R2;
    S = (kron(I,X()^2 * H^2 * X()) * kron(X()^2,I) * (CX * kron(I,R2^8))^3 * kron(X(),I))^2;
    S = removeBorrowedAncilla(S,2);
end

function Sdag = Sdag()
    global I R2;
    Sdag = ( kron(Xdag(),I) * (kron(I,R2)*CXdag())^3 * kron(X(),I) * kron(I,Xdag()*Hdag()^2*X()) )^2;
    Sdag = removeBorrowedAncilla(Sdag,2);
end

% zeta^8 * S; so the dagger of this is zeta S^\dagger
function S8 = S8()
    global H R2;
    S8 = R2 * X()^2 * Hdag()^2 * X() * H^2 * X() * R2^8 * X();
end

% zeta S^\dagger
function Sdag1 = Sdag1()
    global H CX R2;
    Sdag1 = X()^2 * R2 * X()^2 * Hdag()^2 * X()^2 * H^2 * X() * R2^8;
end

% global phase of \zeta = \omega_2 1-qutrit gate
function zgp = zgp()
    zgp = Sdag1() * S();
end

% global phase of \zetadag
function zdaggp = zdaggp()
    zdaggp = Sdag() * S8();
end

% X gate controlled on just one basis state
function ZCX = ZCX(ctrl)
    global I;
    if nargin < 1
        ctrl = 0;
    end
    ZCX = kron(X()^ctrl,I) * ZCXutgp()*kron(zdaggp()^2,I) * kron(Xdag()^ctrl,I);
end

function ZCXdag = ZCXdag(ctrl)
    global d;
    if nargin < 1
        ctrl = 0;
    end
    ZCXdag = ZCX(ctrl)^(d-1);
end

% X_[00,01] level operator, i.e. |0>-controlled X_[0,1] gate, for d = 3
% which is unitary.  For d > 3, known constructions have computational
% basis initialization and deterministic measurement due to X_[0,1].
function X00_01 = X00_01()
    global d I CX;
    udZCX = reIndex(ZCX(),[2 1]); % upside-down ZCX gate
    udZCXdag = reIndex(ZCXdag(),[2 1]);
    perm00_10_01 = udZCX * ZCX() * udZCXdag * ZCXdag();
    udCX = reIndex(CX,[2 1]);
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

% The \omegadag_[22] level operator up to a global phase of \zeta.
% Is the |2>-controlled \zeta Sdag gate.
function wdag1_22 = wdag1_22()
    global I H R2;
    wdag1_22 = ZCX(2)^2 * kron(I,R2*X()^2*Hdag()^2) * ZCX(2)^2 * ...
        kron(I,H^2*X()*R2^8);
end

% The \zeta_[2] level operator
function z2 = z2()
    z2 = removeBorrowedAncilla(wdag1_22()*w22(),2);
end

% The H_[20,21,22] level operator up to controlled global phase -1
function H2m = H2m()
    global I H;
    X0_2 = X()*Hsqm()*X()^2;
    S2_2perm=kron(I,X0_2)*wdag1_22()^2*kron(I,X0_2);
    H2m = S2_2perm*kron(I,Hdag())*S2_2perm*kron(I,H)*S2_2perm;
end

function Hdag2m = Hdag2m()
    global I H;
    X0_2 = X()*Hsqm()*X()^2;
    Sdag1_2perm=kron(I,X0_2)*wdag1_22()*kron(I,X0_2);
    Hdag2m = Sdag1_2perm*kron(I,Hdag())*Sdag1_2perm*kron(I,H)*Sdag1_2perm;
end

% The (-1)_[2] level operator
function m2 = m2()
    global I;
    X21_22 = kron(Xdag(),X())*X00_01()*kron(X(),Xdag());
    m2 = H2m()^2 * X21_22 * kron(Sdag(),I);
    m2 = removeBorrowedAncilla(m2,2);
end

% The H_[20,21,22] 2-qutrit level operator
function H2 = H2()
    global I;
    H2 = H2m() * kron(m2(),I);
end

% The Hdag_[20,21,22] 2-qutrit level operator
function Hdag2 = Hdag2()
    global I;
    Hdag2 = Hdag2m() * kron(m2(),I);
end