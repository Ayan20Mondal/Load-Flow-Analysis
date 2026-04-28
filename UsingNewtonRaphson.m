% | From | To | R | X | B/2 | X'mer |
% | Bus | Bus | | | | |
linedata = [ 1 2 0.10 0.20 0.02 1;
1 4 0.05 0.20 0.02 1;
1 5 0.08 0.30 0.03 1;
2 3 0.05 0.25 0.03 1;
2 4 0.05 0.10 0.01 1;
2 5 0.10 0.30 0.02 1;
2 6 0.07 0.20 0.025 1;
3 5 0.12 0.26 0.025 1;
3 6 0.02 0.10 0.01 1;
4 5 0.20 0.40 0.04 1;
5 6 0.10 0.30 0.03 1];
fb = linedata(:,1); % From bus number...
tb = linedata(:,2); % To bus number...
r = linedata(:,3); % Resistance, R...
x = linedata(:,4); % Reactance, X...
b = linedata(:,5); % Ground Admittance, B/2...
a = linedata(:,6); % Tap setting value..
z = r + i*x; % Z matrix...
y = 1./z; % To get inverse of each element...
b = i*b; % Make B imaginary...
nb = max(max(fb),max(tb)); % no. of buses...
nl = length(fb); % no. of branches...
Y = zeros(nb,nb); % Initialise YBus...
% Formation of the Off Diagonal Elements...
for k = 1:nl
Y(fb(k),tb(k)) = Y(fb(k),tb(k)) - y(k)/a(k);
Y(tb(k),fb(k)) = Y(fb(k),tb(k));
end
% Formation of Diagonal Elements....
for m = 1:nb
for n = 1:nl
if fb(n) == m
Y(m,m) = Y(m,m) + y(n)/(a(n)^2) + b(n);
elseif tb(n) == m
Y(m,m) = Y(m,m) + y(n) + b(n);
end
end
end
Y
% |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi | Qmin | Qmax |
busd = [ 1 1 1.05 0 0.0 0 0 0 0 0;
2 2 1.05 0 0.5 0 0 0 -0.5 1.0;
3 2 1.07 0 0.6 0 0 0 -0.5 1.5;
4 3 1.0 0 0.0 0 0.7 0.7 0 0;
5 3 1.0 0 0.0 0 0.7 0.7 0 0;
6 3 1.0 0 0.0 0 0.7 0.7 0 0];
BMva = 100; % Base MVA..
bus = busd(:,1); % Bus Number..
type = busd(:,2); % Type of Bus 1-Slack, 2-PV, 3-PQ..
V = busd(:,3); % Specified Voltage..
del = zeros(length(V),1); % Voltage Angle..
Pg = busd(:,5)/BMva; % PGi..
Qg = busd(:,6)/BMva; % QGi..
Pl = busd(:,7)/BMva; % PLi..
Ql = busd(:,8)/BMva; % QLi..
Qmin = busd(:,9)/BMva; % Minimum Reactive Power Limit..
Qmax = busd(:,10)/BMva; % Maximum Reactive Power Limit..
nbus = max(bus); % To get no. of buses..
P = Pg - Pl; % Pi = PGi - PLi..
Q = Qg - Ql; % Qi = QGi - QLi..
Psp = P; % P Specified
Qsp = Q; % Q Specified
G = real(Y); % Conductance..
B = imag(Y); % Susceptance..
pv = find(type == 2 | type == 1); % Index of PV Buses..
pq = find(type == 3); % Index of PQ Buses..
npv = length(pv); % Number of PV buses..
npq = length(pq); % Number of PQ buses..
Tol = 1;
Iter = 1; % iteration starting
while (Tol > 1e-5) % Iteration starting..
P = zeros(nbus,1);
Q = zeros(nbus,1);
% Calculate P and Q
for i = 1:nbus
for k = 1:nbus
P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
end
end
% Checking Q-limit violations..
if Iter <= 7 && Iter > 2 % Only checked up to 7th iterations..
for n = 2:nbus
if type(n) == 2
QG = Q(n)+Ql(n);
if QG < Qmin(n)
V(n) = V(n) + 0.01;
elseif QG > Qmax(n)
V(n) = V(n) - 0.01;
end
end
end
end
% Calculate change from specified value
dPa = Psp-P;
dQa = Qsp-Q;
k = 1;
dQ = zeros(npq,1);
for i = 1:nbus
if type(i) == 3
dQ(k,1) = dQa(i);
k = k+1;
end
end
dP = dPa(2:nbus);
M = [dP; dQ]; % Mismatch Vector
% Jacobian
% J1 - Derivative of Real Power Injections with Angles..
J1 = zeros(nbus-1,nbus-1);
for i = 1:(nbus-1)
m = i+1;
for k = 1:(nbus-1)
n = k+1;
if n == m
for n = 1:nbus
J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
end
J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
else
J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
end
end
end
% J2 - Derivative of Real Power Injections with V..
J2 = zeros(nbus-1,npq);
for i = 1:(nbus-1)
m = i+1;
for k = 1:npq
n = pq(k);
if n == m
for n = 1:nbus
J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
end
J2(i,k) = J2(i,k) + V(m)*G(m,m);
else
J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
end
end
end
% J3 - Derivative of Reactive Power Injections with Angles..
J3 = zeros(npq,nbus-1);
for i = 1:npq
m = pq(i);
for k = 1:(nbus-1)
n = k+1;
if n == m
for n = 1:nbus
J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
end
J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
else
J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
end
end
end
% J4 - Derivative of Reactive Power Injections with V..
J4 = zeros(npq,npq);
for i = 1:npq
m = pq(i);
for k = 1:npq
n = pq(k);
if n == m
for n = 1:nbus
J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
end
J4(i,k) = J4(i,k) - V(m)*B(m,m);
else
J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
end
end
end
J = [J1 J2; J3 J4]; % Jacobian Matrix
X = J\M; % INV(J) x M, Correction Vector..
dTh = X(1:nbus-1); % Change in Voltage Angle..
dV = X(nbus:end); % Change in Voltage Magnitude..
% Update State Vectors (Voltage Angle & Magnitude)
del(2:nbus) = dTh + del(2:nbus);
k = 1;
for i = 2:nbus
if type(i) == 3
V(i) = dV(k) + V(i);
k = k+1;
end
end
Iter = Iter + 1;
Tol = max(abs(M));
end
Del = 180/pi*del; % Convert radians to degrees
nb = length(V); % Number of buses
Iij = zeros(nb,nb);
Sij = zeros(nb,nb);
% Bus Current Injections..
I = Y*V;
Im = abs(I);
Ia = angle(I);
%Line Current Flows..
for m = 1:nb
for n = 1:nl
if fb(n) == m
p = tb(n);
Iij(m,p) = -(V(m) - V(p)*a(n))*Y(m,p)/a(n)^2 + b(n)/a(n)^2*V(m);
Iij(p,m) = -(V(p) - V(m)/a(n))*Y(p,m) + b(n)*V(p);
elseif tb(n) == m
p = fb(n);
Iij(m,p) = -(V(m) - V(p)/a(n))*Y(p,m) + b(n)*V(m);
Iij(p,m) = -(V(p) - V(m))*Y(m,p)/a(n)^2 + b(n)/a(n)^2*V(p);
end
end
end
Iijr = real(Iij);
Iiji = imag(Iij);
% Line Power Flows..
for m = 1:nb
for n = 1:nb
if m ~= n
Sij(m,n) = V(m)*conj(Iij(m,n))*BMva;
end
end
end
Pij = real(Sij);
Qij = imag(Sij);
% Line Losses..
Lij = zeros(nl,1);
for m = 1:nl
p = fb(m); q = tb(m);
Lij(m) = Sij(p,q) + Sij(q,p);
end
Lpij = real(Lij);
Lqij = imag(Lij);
% Bus Power Injections..
Si = zeros(nb,1);
for i = 1:nb
for k = 1:nb
Si(i) = Si(i) + conj(V(i))* V(k)*Y(i,k)*BMva;
end
end
Pi = real(Si)
Qi = -imag(Si)
Pg = Pi+Pl
Qg = Qi+Ql
Pl
Ql
V
Del