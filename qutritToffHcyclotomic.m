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

% testU(X);
% disp(H);
% testU(CCX);
% testU(Xdag());
% disp(Hdag());
% testU(H2m());
% testU(CX());
% testU(CXdag());
% testU(Swapm());
% 
% % The constructions from here on are specific to d = 3. For d > 3 it's
% % easier to use the gate set {X, H, ZCX} instead of {X, H, CCX}.
% testU(ZCX());
% testU(ZCXdag());
% testU(X00_01());
% testU(X000_001());
% testU(S());
testU(w22());
1;

% P3 = swap*ZCX*swap * ZCX * swap*ZCX^(d-1)*swap * ZCX^(d-1);
% Z_OCX01 = P3*swap*CX*swap*P3*swap*CX^(d-1)*swap;
% X01 = kron(k0',I)*Z_OCX01*kron(k0,I);
% ZCX01 = kron(X,X01) * (Z_OCX01* kron(X^2,I))^((d-1)/2);
% CX01 = kron(X,I) * (ZCX01*kron(X^2,I))^((d-1)/2);
% ZZCX01 = kron(I,CX01) * kron(ZCX,I) * kron(I,CX01) * kron(ZCX^(d-1),I) * kron(swap,I)*kron(I,ZCX01)*kron(swap,I);
% ZZCX = (ZZCX01 * kron(II,X))^(d-1) * kron(II,X);
% 
% Q0 = kron(I,k0'*X^(d-1)*H) * ZCX * kron(I,H^3*X*k0);

% %if d = 3
%     ZZCw = (kron(II,H^3)*ZZCX*kron(II,H) * kron(II,X^(d-1) * X01 * X))^2;
%     %ZZCwdag = (kron(II,H)*ZZCX*kron(II,H^3) * kron(II,X^(d-1) * X01 * X))^2;
%     ZCS = removeBorrowedAncilla(kron(kron(I,X^(d-1)),I) * ZZCw * kron(kron(I,X),I),1);
%     TCiHdag = kron(X^(d-1),H^2) * ZCS * kron(I,H^2) * ZCS * kron(I,H^3) * ZCS * kron(I,H^2) * ZCS * kron(I,H) * ZCS * kron(I,H^2) * ZCS * kron(X,I);
%     %TCminusiH = kron(X^(d-1),H^2) * ZCS^2 * kron(I,H^2) * ZCS^2 * kron(I,H^3) * ZCS^2 * kron(I,H^2) * ZCS^2 * kron(I,H) * ZCS^2 * kron(I,H^2) * ZCS^2 * kron(X,I);
%     RtensI = TCiHdag^2 * kron(X^(d-1),X)*ZCX01*kron(X,X^(d-1));
%     R = removeBorrowedAncilla(RtensI,1);
% %end

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

function Xdag = Xdag()
    global d X;
    Xdag = X^(d-1);
end

% H^\dagger with exact global phase. Is H^3 up to a global phase.
function Hdag = Hdag()
    global d H;
    Hdag = H^(4*d-1);
end

% Sends |x> to -|-x mod d>. Is H^2 up to a global phase.
function H2m = H2m()
    global d H;
    H2m = H^(2*d);
end

function CX = CX()
    global d II X H CCX;
    CX = (kron(X * H2m() * Xdag(), II) * CCX^((d+1)/2))^2;
    CX = removeBorrowedAncilla(CX,1);
end

function CXdag = CXdag()
    global d;
    CXdag = CX()^(d-1);
end

% Swap gate, with a global phase of -1.
function Swapm = Swapm()
    global I H;
    CXupsidedown = reIndex(CX(), [2 1]);
    Swapm = kron(I,H^2)*CX * kron(H^2,I)*CXupsidedown * kron(I,H^2)*CX;
end

% |0>-controlled X gate: Does on target X if control is |0>, else does I.
function ZCX = ZCX()
    global d I X H CCX;
    nZCX = kron(kron(I,Hdag()^2),I) * CCX * kron(kron(I,H^2),I) * kron( ...
        CXdag(),I) * CCX * kron(CX(),I); % |not 0>-controlled X
    nZCX = removeBorrowedAncilla(nZCX,2);
    ZCX = kron(I,X) * nZCX^(d-1);
end

% |0>-controlled Xdag: Does on target Xdag if control is |0>, else does I.
function ZCXdag = ZCXdag()
    global d;
    ZCXdag = ZCX()^(d-1);
end

% X_[00,01] level operator, i.e. |0>-controlled X_[0,1] gate. Global phase
% -1. This is for d = 3 which is unitary.  For d > 3, known constructions
% have computational basis initialization and deterministic measurement
% to do X_[0,1], unless that can be done unitarily in the gate set.
function X00_01 = X00_01()
    global d I X H;
    udZCX = reIndex(ZCX(),[2 1]);
    udZCXdag = reIndex(ZCXdag(),[2 1]);
    perm_00_10_01 = udZCX * ZCX() * udZCXdag * ZCXdag();
    udCX = reIndex(CX(),[2 1]);
    udCXdag = reIndex(CXdag(),[2 1]);
    X_0_1 = X^2*H2m()*X; % This holds only for d = 3
    X00_01 = kron(X,X_0_1) * (perm_00_10_01*udCX*perm_00_10_01* ...
        udCXdag*kron(X^2,I))^((d-1)/2);
end

function X000_001 = X000_001()
    global d I X;
    CX0_1 = (kron(Xdag()^2,I) * X00_01)^((d-1)/2) * kron(Xdag(),I);
    X000_001 = kron(I,CX0_1) * kron(ZCX(),I) * kron(I,CX0_1) * ...
        kron(ZCXdag(),I) * reIndex(X00_01(),[1 3 2]);
end

% The 1-qutrit S gate. Is the \omega_[2] operator applying \omega to |2>.
function S = S()
    global X H;
    S = (kron(X^2,H)*ZCX()*kron(X,H^3*X^2*H^2*X))^2;
    S = removeBorrowedAncilla(S,2);
end

% The \omega_[22] level operator which applies phase \omega to |22>.
% Is the |2>-controlled S gate.
function w22 = w22()
    global d II X H;
    ZZCX = kron(II,X)*(kron(II,X)*X000_001())^(d-1); % |00>-controlled X, i.e. X_[000,001,002]
    w22 = (kron(kron(X^2,X^2),H)*ZZCX*kron(kron(X,X),H^3*X^2*H^2*X))^2;
    w22 = removeBorrowedAncilla(w22,3);
end
