function res = verifyMatrixContents(A, B, tol)
%VERIFYMATRIXCONTENTS Check whether two matrices contain the same data.
%   Since floating point operations are imprecise, a tolerance is allowed.

% If the dimensions aren't the same, don't even look at the contents.
if (~isequal(size(A), size(B)))
    res = 0;
    return;
end

capacity = numel(A);

res = sum(abs(A(:) - B(:)) <= tol) == capacity;

end
