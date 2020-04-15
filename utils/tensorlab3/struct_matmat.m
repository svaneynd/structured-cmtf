function [x,state] = struct_matmat(z,task)
%STRUCT_MATMAT Compute matrix matrix product. 
%   [X,STATE] = STRUCT_MATMAT(Z,TAKS]) computes the matrix matrix product
%   Z{1}*Z{2}*...*Z{N} in which N = length(Z).  The structure state stores
%   information which is reused in computing the right and left Jacobian-vector
%   products.
%
%   STRUCT_MATMAT(Z,T) computes the right or left Jacobian-vector product of
%   this transformation, depending on the structure task. Use the structure
%   state and add the field 'r' of the same shape as z or the field 'l' of the
%   same shape as x to obtain the structure task for computing the right and
%   left Jacobian-vector products
%   
%      (dF(:)/dz(:).')*task.r(:) and
%      (dF(:)/dz(:).')'*task.l(:) + conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%   
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial
%   derivative which treats conj(z) (z) as constant. The output has the
%   same shape as x or z for the right and left Jacobian-vector products,
%   respectively.

%   Authors: Nico Vervliet       (Nico.Vervliet@kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

    if nargin < 2, task = []; end

    state = [];
    left   = isstruct(task) && isfield(task, 'l') && ~isempty(task.l);
    right  = isstruct(task) && isfield(task, 'r') && ~isempty(task.r);
    
    if ~left && ~right
        if ~iscell(z)
            x = z;
        else
            x = z{1};
            for n = 2:numel(z)
                x = x * z{n};
            end 
        end 
    elseif left
        x = cell(size(z));
        for n = 1:numel(z)
            tmp = 1;
            for m = 1:n-1
                tmp = z{m}' * tmp;
            end
            tmp = tmp * task.l;
            for m = numel(z):-1:n+1
                tmp = tmp * z{m}';
            end 
            x{n} = tmp;
        end 
    elseif right
        x = 0;
        for n = 1:numel(z)
            tmp = 1;
            for m = 1:numel(z)
                if m == n
                    tmp = tmp * task.r{n};
                elseif m < n 
                    tmp = tmp * z{m};
                else 
                    tmp = tmp*z{m};
                end 
            end 
            x = x + tmp;
        end 
    end 
end
