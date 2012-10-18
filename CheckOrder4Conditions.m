function [methodorder, conditions] = CheckOrder4Conditions(a,b)
% Jonathan Ong
% jonathanong@berkeley.edu
%
% This function checks the order conditions of a Runge-Kutta method up to
% its order 4 conditions with inputs A and b from the method's butcher
% array.
%
% Outputs:
% methodorder - tells you the method's order up to order 4 (0,1,2,3,4)
% conditions - 3 columns matrix:
%               [order conditions, calculated order conditions, difference]
%
% Ex. Explicit Runge-Kutta of Order 4
% a = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
% b = [1 2 2 1]/6;
% Ex. Hammer-Holllingsworth of Order 4
% a = [1/4 1/4-sqrt(3)/6; 1/4+sqrt(3)/6 1/4];
% b = [1 1]/2;
% Ex. Butcher's Lobatto formula of order 4
% a = [0 0 0; 1/4 1/4 0; 0 1 0];
% b = [1/6 2/3 1/6];

% Order(1) = Order 1 conditions = 1
% Order(2) = Order 2 conditions = 1/2
% Order(3) = Order 3a conditions = 1/3
% Order(4) = Order 3b conditions = 1/6
% Order(5) = Order 4a conditions = 1/8
% Order(6) = Order 4b conditions = 1/4
% Order(7) = Order 4c conditions = 1/24
% Order(8) = Order 4d conditions = 1/12

if nargin < 2
    error('Insufficient inputs')
end

[n,m] = size(a);
c = length(b);
if n ~= m
    error('A of butcher array must be square')
end
if n ~= c
    error('Dimensions of A of butcher array must be equal to the dimension of b')
end

conditions = [1;1/2;1/3;1/6;1/8;1/4;1/24;1/12];
Order = zeros(8,1);

% Order 1
for l = 1:c
    Order(1) = Order(1) + b(l);
end

% Order 2
for l = 1:c; for i = 1:c;
        Order(2) = Order(2) + b(l)*a(l,i);
    end; end;

% Order 3
for l = 1:c; for i = 1:c; for j = 1:c;
            Order(3) = Order(3) + b(l)*a(l,i)*a(l,j);
        end; end; end;

% Order 4
for l = 1:c; for i = 1:c; for j = 1:c;
            Order(4) = Order(4) + b(l)*a(l,i)*a(i,j);
        end; end; end;

% Order 5
for l = 1:c; for i = 1:c; for j = 1:c; for k = 1:c;
                Order(5) = Order(5) + b(l)*a(l,i)*a(l,j)*a(j,k);
            end; end; end; end;

% Order 6
for l = 1:c; for i = 1:c; for j = 1:c; for k = 1:c;
                Order(6) = Order(6) + b(l)*a(l,i)*a(l,j)*a(l,k);
            end; end; end; end;

% Order 7
for l = 1:c; for i = 1:c; for j = 1:c; for k = 1:c;
                Order(7) = Order(7) + b(l)*a(l,i)*a(i,j)*a(j,k);
            end; end; end; end;

% Order 8
for l = 1:c; for i = 1:c; for j = 1:c; for k = 1:c;
                Order(8) = Order(8) + b(l)*a(l,i)*a(i,j)*a(i,k);
            end; end; end; end;

Difference = Order - conditions;
conditions = [Order, conditions, Difference];

NotSat = find(abs(Difference) > 1e-15,1,'first');
if NotSat == 1
    methodorder = 0;
elseif NotSat == 2
    methodorder = 1;
elseif NotSat == 3 || NotSat == 4
    methodorder = 2;
elseif NotSat == 5 || NotSat == 6 || NotSat == 7 || NotSat == 8
    methodorder = 3;
elseif isempty(NotSat)
    methodorder = 4;
else
    methodorder = -1;
end