function [c, ceq] = totalTimeConstraint(x)
c = [sum(x([1 4 ])) + 4 .* 7.000000 - 57.000000];
ceq = [];
end