function [c, ceq] = totalTimeConstraint(x)
c = [sum(x([1 2 7 8 ])) + 4 .* 7.000000 - 96.000000];
ceq = [];
end