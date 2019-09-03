function [x, y, z] = readFile(fileName)
  A = dlmread(fileName, '\t');
  
  x_tmp = A(1, 1:end - 1);
  y_tmp = A(2:end, 1);
  
  A(1, :) = [];
  A(:, 1) = [];
  
  [r, c] = size(A);
  
  x = zeros(r, c);
  y = zeros(r, c);
  z = zeros(r, c);
  for i = 1 : r
    x(i, :) = x_tmp(:);
  end
  for j = 1 : c
    y(:, j) = y_tmp(:);
  end
  z = A;
%  for i = 1 : r
%    for j = 1 : c
%      x(i, :) = x_tmp(:);
%      y(i, j) = y_tmp(i);
%      z(i, j) = A(i, j);
%    end
%  end
end