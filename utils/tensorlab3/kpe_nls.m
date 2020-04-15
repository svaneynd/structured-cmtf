function [x,output] = kpe_nls(A, b, x0, varargin)
%KPE_NLS Kronecker product equations by nonlinear least squares. 
%   [x,output] = kpe_nls(A, b, x0) computes the solution of the system Ax=b
%   where x has a Kronecker product structure, i.e., x = kr(x{N}, x{N-1}, ...,
%   x{1}). The matrix A can be given as a matrix of size I_0 x I_1I_2...I_N, or
%   as a tensor of size I_0 x I_1 x I_2 x ... x I_N. The algorithm is
%   initialized by the Kronecker product factors in x0. The structure output
%   returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   kpe_nls(A, b, x0, options) may be used to set the following options:
%
%      options.Algorithm =   - The desired optimization method.
%      [@nls_gncgs| ...
%       {@nls_gndl}|@nls_lm]
%      options.M =           - The preconditioner to use when
%      [{'block-Jacobi'}|...   options.LargeScale is true.
%       false]
%      options.<...>         - Parameters passed to the selected method,
%                              e.g., options.TolFun, options.TolX and
%                              options.PlaneSearchOptions. See also help
%                              [options.Algorithm].

%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Convert A to a tensor
    N = length(x0);
    if ndims(A) ~= N + 1
        if prod(cellfun(@length, x0)) ~= size(A, 2), 
            error('kpe_nls:dimensionMismatch', 'Dimensions don''t agree');
        end
        A = reshape(A, [size(A,1); cellfun(@length, x0(:))]');
    end
    nbrows = size(A,1);

    % permute the first mode to the middle to minimize middle contractions
    m1 = ceil(0.5*(N+1));
    A = permute(A, [2:m1 1 m1+1:N+1]);
    modes = [1:m1-1 m1+1:N+1];
    
    % Create helper functions and variables
    isfunc = @(f)isa(f,'function_handle');
    xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
    funcs = {@nls_gndl,@nls_gncgs,@nls_lm};
    
    % Parse options
    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('Algorithm', funcs{find(cellfun(xsfunc,funcs),1)});
    p.addParameter('CGMaxIter', 10);
    p.addParameter('Display', 0);
    p.addParameter('M', 'block-Jacobi');
    p.parse(varargin{:});
    options = [fieldnames(p.Results)'  fieldnames(p.Unmatched)';
               struct2cell(p.Results)' struct2cell(p.Unmatched)'];
    options = struct(options{:});
        
    % Cache variables
    cache.offset = cellfun(@numel,x0(:));
    cache.offset = cumsum([1; cache.offset]);
    
    % Call the optimization method.
    dF.JHJx = @JHJx;
    dF.JHF = @grad;
    switch options.M
      case 'block-Jacobi', dF.M = @M_blockJacobi;
      otherwise, if isfunc(options.M), dF.M = options.M; end
    end
    [x,output] = options.Algorithm(@objfun,dF,x0(:).',options);
    output.Name = func2str(options.Algorithm);

    function fval = objfun(z)
    % KPE objective function
        fval = contract(A, z, modes) - b;
        cache.residual = fval;
        fval = 0.5*(fval(:)'*fval(:));
    end

    function grad = grad(z)
    % KPE scaled conjugate cogradient.
        E = cache.residual;
        offset = cache.offset;
        grad = nan(offset(end)-1,1);
        
        for n = 1:N
            cache.J{n} = contract(A, z([1:n-1 n+1:N]), modes([1:n-1 n+1:N]));
            if n > N/2 % because of permuting first mode to middle
                cache.J{n} = cache.J{n}.';
            end
            grad(offset(n):offset(n+1)-1) = conj(cache.J{n})*E;
            cache.JHJinv{n} = pinv(conj(cache.J{n})*cache.J{n}.');
        end
    end

    function y = JHJx(~,y)
    % KPE fast Gramian vector product
        offset = cache.offset;
        J  = cache.J;
        Jx = zeros(nbrows, 1);
        for n = 1:N
            Jx = Jx + (y(offset(n):offset(n+1)-1).'*J{n}).';
        end
        for n = 1:N
            y(offset(n):offset(n+1)-1) = conj(J{n}) * Jx;
        end
    end

    function y = M_blockJacobi(~,q)
    % Solve M*y = q, where M is a block-diagonal approximation for JHJ.
        y = nan(size(q));
        offset = cache.offset;
        for n = 1:N
            idx = offset(n):offset(n+1)-1;
            y(idx) = cache.JHJinv{n}*q(idx);
        end
    end
   
end