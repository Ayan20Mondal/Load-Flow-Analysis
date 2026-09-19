clc;
clear;
%% ------------------ INPUT DATA ------------------
nb = 3; % number of buses
% Bus data: [Bus No, Type(1=Slack,2=PV,3=PQ), V, delta, P, Q]
busdata = [1 1 1.06 0 0 0; % Slack
2 3 1.00 0 -0.5 -0.2; % PQ
3 3 1.00 0 -0.6 -0.25; ]; % PQ
% Line data: [From To R X]
linedata = [1 2 0.02 0.06;
1 3 0.08 0.24;
2 3 0.06 0.18;];
%% ------------------ Y-BUS FORMATION ------------------
Y = zeros(nb);
for k = 1:size(linedata,1)
i = linedata(k,1);
j = linedata(k,2);
z = linedata(k,3) + 1i*linedata(k,4);
y = 1/z;
Y(i,i) = Y(i,i) + y;
Y(j,j) = Y(j,j) + y;
Y(i,j) = Y(i,j) - y;
Y(j,i) = Y(i,j);
end
%% ------------------ INITIALIZATION ------------------
V = busdata(:,3) .* exp(1i*deg2rad(busdata(:,4)));
P = busdata(:,5);
Q = busdata(:,6);
tol = 1e-6;
max_iter = 100;
iter = 0;
%% ------------------ GAUSS-SEIDEL ITERATION ------------------
while iter < max_iter
V_old = V;
for i = 1:nb
type = busdata(i,2);
if type == 1
continue; % Slack bus
end
sumYV = 0;
for j = 1:nb
if j ~= i
sumYV = sumYV + Y(i,j)*V(j);
end
end
% PQ bus update
V(i) = (1/Y(i,i)) * ((P(i) - 1i*Q(i))/conj(V(i)) - sumYV);
end
% Convergence check
if max(abs(V - V_old)) < tol
break;
end
iter = iter + 1;
end
%% ------------------ RESULTS ------------------
disp('Bus Voltages (p.u.):');
for i = 1:nb
fprintf('Bus %d: |V|=%.4f, Angle=%.2f deg\n', ...
i, abs(V(i)), rad2deg(angle(V(i))));
end
fprintf('\nTotal Iterations = %d\n', iter);