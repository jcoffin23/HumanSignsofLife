function [outty] = threedto2d(mat3)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
  A = [size(mat3)];
  A = A(1:2);
outty = zeros(size(mat3,3),A(1) * A(2));
for k = 1: size(mat3,2) * size(mat3,1)

  [I,J] = ind2sub(A ,k);
  outty(:,k) = mat3(I,J,:);
end

end

