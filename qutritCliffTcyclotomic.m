%% Implementation of qutrit controlled Clifford+T gates
% For the paper Exact synthesis of multiqutrit Clifford-cyclotomic circuits
% An earlier implementation, with fewer gates and where some of the qutrit
% Clifford gates differ in definition by a global phase, is at https://git
% hub.com/lia-approves/qudit-circuits/tree/main/qutrit_control_Clifford_T

%% declare and define global variables of the minimal gate set
% global I II w s h cx x01 x02 x12 swap t;
global CX H R2 I II;

def_gates();
1;

function def_gates()
    global CX H R2 I II;
    w = exp(i*2*pi/3);
    z = exp(i*2*pi/9);
    I = eye(3);     % single-qutrit identity gate, i.e. 3 x 3 identity matrix
    II = eye(3^2);  % two-qutrit identity gate
    R2 = [1 0 0; 0 z 0; 0 0 z^2];  % R2 // Def.?
    H = -w^2/sqrt(-3)*[1 1 1; 1 w w^2; 1 w^2 w^4];   % H // Def.?
    zero = zeros([3]);  % 3 x 3 zero matrix
    CX = [I zero zero; zero X() zero; zero zero X()^2];   % CX, i.e. tau(10 11 12)(20 22 21) // Def.?
end

function Z = Z()
    global R2;
    Z = R2^3;
end

function X = X()
    global H;
    X = H^7 * Z() * H^5;
end

function S = S()
    global CX H R2 I;
    S = (kron(I,X()^2 * H^2 * X()) * kron(X()^2,I) * (CX * kron(I,R2^8))^3 * kron(X(),I))^2;
    S = removeBorrowedAncilla(S);
end

% negative X_[a,b]
function NXab = NXab(a,b)
    global H;
    if any(mod(a-b,3) == [1,2])
        NXab = X()^(a+b)^2 * H^2 * X()^(a+b);
    else
        ME = MException('component:invalidAB',sprintf('Invalid a=%d, b=%d',a,b));
        throw(ME)
    end
end

% X gate controlled on just one basis state
function BCX = BCX(ctrl)
    global CX H R2 I;
    BCX = kron(X^ctrl,H^11) * (kron(I,R2) * CX)^3 * kron(NXab(0,2) * S() * NXab(0,2) * X^2^ctrl,H);
end

% zeta^8 * S; so the dagger of this is zeta S^\dagger
function S8 = S8()
    global CX H R2;
    S8 = H^8 * R2 * X()^2 * H^2 * X() * H^2 * X() * R2^8 * X();
end

% zeta S^\dagger
function Sdag1 = Sdag1()
    global CX H R2;
    Sdag1 = X()^2 * R2 * X()^2 * H^2 * X()^2 * H^2 * X * R2^8 * H^8;
end

function V = removeBorrowedAncilla(U)
    global I;
    V = U(1:3:end,1:3:end);
    if any(any(abs(kron(V,I) - U) > 1e-5))
        str = sprintf([repmat('%.3f ', 1, size(U, 2)) '\n'], U);
        ME = MException('component:noBorrowedAncilla',[str 'has no borrowed ancilla'],str);
        throw(ME)
    end
end