function [c, ceq] = totalTimeConstraint(x)
c = [sum(x([])) + 9 .* 7.000000 - 96.000000];
ceq = [];
end