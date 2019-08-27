clc;
clear all;
close all;

D = dlmread('det.csv', '\t');
TR = dlmread('trace.csv', '\t');

omega = D(1, 1:end-1);
e = D(2:end, 1);

D(1, :) = [];
D(:, 1) = [];
TR(1, :) = [];
TR(:, 1) = [];

fprintf(1, "D in [%.15e, %.15e]\n", min(min(D)), max(max(D)));

[ee, oo] = meshgrid(e, omega);

%contour(oo, ee, TR, [2 2]);

contour(oo, ee, TR, [2 2]);

eps = 1;
omega = 0.31;
eqs = @(t, q)[q(2) ; -omega^2 * (1 + eps * cos(t)) * q(1)];
opts = odeset('MaxStep', 1e-2, 'AbsTol', 1e-6, 'RelTol', 1e-6, 'InitialStep', 1e-2);
[t, q] = ode45(eqs, [0, 50*pi], [1e-4, 0], opts);

f2 = figure(2);
subplot(2, 1, 1);
plot(t, q(:, 1));
grid on
subplot(2, 1, 2);
plot(t, q(:, 2));
grid on
