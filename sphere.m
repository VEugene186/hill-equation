clc;
clear all;
close all;

A = dlmread('solution.csv', '\t');


g = norm(A(1, 2:4));

phi = linspace(0, 2*pi, 41);
theta = linspace(0, pi, 21);
[pp, tt] = meshgrid(phi, theta);
x = 0.999 * g * cos(pp) .* sin(tt);
y = 0.999 * g * sin(pp) .* sin(tt);
z = 0.999 * g * cos(tt);

mesh(x, y, z);
hold on
plot3(A(:, 2), A(:, 3), A(:, 4), 'ko', 'markersize', 1);
axis('square');
