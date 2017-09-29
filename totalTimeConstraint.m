function [c, ceq] = totalTimeConstraint(x)
c = [sum(x([1 2 3 4 ])) + 5 .* 7.000000 - 96.000000];
ceq = [];
end