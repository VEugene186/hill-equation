clc;
clear all;
close all;

[omega, J, tr] = readFile('trace.csv');
f1 = figure(1);
h = semilogx(1, 1);
hold on;
contour(omega, J, tr, [2 2], 'linecolor', 'k');
xlabel('Omega');
ylabel('J');
grid on
delete(h);

return

D = dlmread('det.csv', '\t');
TR = dlmread('trace.csv', '\t');

omega = D(1, 1:end-1);
e = D(2:end, 1);

D(1, :) = [];
D(:, 1) = [];
TR(1, :) = [];
TR(:, 1) = [];

fprintf(1, "D in [%.15e, %.15e]\n", min(min(D)), max(max(D)));

[oo, ee] = meshgrid(omega, e);


D1 = dlmread('det1.csv', '\t');
TR1 = dlmread('trace1.csv', '\t');

omega1 = D1(1, 1:end-1);
e1 = D1(2:end, 1);

D1(1, :) = [];
D1(:, 1) = [];
TR1(1, :) = [];
TR1(:, 1) = [];

fprintf(1, "D in [%.15e, %.15e]\n", min(min(D)), max(max(D)));

[oo, ee] = meshgrid(omega, e);
[oo1, ee1] = meshgrid(omega1, e1);
%ee = ee.^2;

f1 = figure(1);
%contour(oo, ee, TR, [2 2]);
h = semilogx(1, 1);
hold on
%set(gca, 'xscale', 'log')
%contour(oo, ee, TR, [2 2], 'linewidth', 2, 'linecolor', 'k');
contour(oo1, ee1, TR1, [2 2], 'linecolor', 'k');
hold on;
contour(oo, ee, TR, [2 2], 'linecolor', 'k');
xlabel('Omega');
ylabel('J');
grid on
delete(h);

text(0.22, 2.5, 'U');
text(0.07, 1.5, 'U');
text(0.04, 1.2, 'U');
text(0.2, 0.8, 'S');
text(0.1, 0.8, 'S');

I1 = 2;
I2 = 3;
I3 = 4;
K = sqrt((I1 - I3) * (I1 - I2) / I1^2 / I2 / I3);
plot([K , K / 2 , K / 3, K / 4], [0 0 0 0], 'ko', 'markerfacecolor', 'k');

FS = findall(f1, '-property', 'FontSize');
set(FS, 'FontSize', 16);
set(FS, 'FontName', 'Times');


%surfc(oo, ee, TR, 'edgecolor', 'none');
%zlim([0 10]);
return
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
