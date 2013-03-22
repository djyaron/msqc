function res = verifyMatrixContents(A, B, tol)
%VERIFYMATRIXCONTENTS Check whether two matrices contain the same data.
%   Since floating point operations are imprecise, a tolerance is allowed.

res = 1;

% If the dimensions aren't the same, don't even look at the contents.
if (~isequal(size(A), size(B)))
    res = 0;
    return;
end

capacity = numel(A);

% Look at each element and compare the absolute difference to the
% tolerance.
for i = 1:capacity
    if (abs(A(i)-B(i)) > tol)
        res = 0;
        return;
    end
end

end

