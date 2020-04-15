function s = cossimil( a , b )
%COSSIMIL computes the cosine similarity between arrays a (of dimension
% NxP) and array b (of dimensions NxQ) and returns this into a PxQ matrix

    if nargin < 2 
        b = a;
    end
    assert(size(a,1)==size(b,1))
    an = bsxfun(@rdivide,a,sqrt(dot(a,a,1)));
    bn = bsxfun(@rdivide,b,sqrt(dot(b,b,1)));
    s = an'*bn;
end

