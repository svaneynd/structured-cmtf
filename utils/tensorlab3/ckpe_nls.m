function [x,output] = ckpe_nls(A, B, x0, varargin)
%CKPE_NLS Coupled Kronecker product equations by nonlinear least squares. 
%  ToDo

% Convert A to a tensor
    N = length(x0);
    if ndims(A) ~= N + 1
%         if prod(cellfun(@length, x0)) ~= size(A, 2), 
        if prod(cellfun(@(x)size(x,1), x0)) ~= size(A, 2), % adapted by SVE because 
            % there was a problem when too strong truncation in some modes
            % was performed:
            % this check would not be correct if the truncated dimension
            % was smaller than the degree of coupling
            error('kpe_nls:dimensionMismatch', 'Dimensions don''t agree');
        end
        A = reshape(A, [size(A,1); cellfun(@length, x0(:))]');
    end
    M = size(B,1);
    Q = size(B,2);
    I = size(x0{2},1);
    nbrows = Q*size(A,1);
    
%     x0 = cellfun(@transpose,x0,'UniformOutput',false);
    
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
%     x = cellfun(@transpose,x,'UniformOutput',false);
    output.Name = func2str(options.Algorithm);

    function fval = objfun(z)
    % KPE objective function
        fval = tens2mat(tmprod(A, {z{1}.',z{2}.'}, [2 3]),1) - B;
        cache.residual = fval(:);
        fval = 0.5*(fval(:)'*fval(:));
    end

    function grad = grad(z)
    % KPE scaled conjugate cogradient.
        E = cache.residual;
        offset = cache.offset;
        grad = nan(offset(end)-1,1);
        
        cache.J{1} = tens2mat(tmprod(A,z{2}.',3),2).'; 
        cache.J{2} = tens2mat(contract(A,z(1),2),1);
        
        grad(offset(1):offset(2)-1) = cache.J{1}'*E;
        cache.JHJinv{1} = pinv(cache.J{1}'*cache.J{1});
        
        grad(offset(2):offset(3)-1) = tens2vec(cache.J{2}'*reshape(E,[M Q]));      
        cache.JHJinv{2} = pinv(cache.J{2}'*cache.J{2}); 
    end

    function y = JHJx(~,y)
    % KPE fast Gramian vector product
        offset = cache.offset;
        J  = cache.J;

        Jx = zeros(nbrows, 1);        
        for n = 1:N
            if n == 1
                Jx = Jx + J{n}*(y(offset(n):offset(n+1)-1));
            else
                % ToDo: do more efficiently
                Jx = Jx + kron(eye(Q),J{n})*(y(offset(n):offset(n+1)-1));
            end
        end
        for n = 1:N
            if n == 1
                y(offset(n):offset(n+1)-1) = J{n}' * Jx;
            else
                % ToDo: do more efficiently
                y(offset(n):offset(n+1)-1) = kron(eye(Q),J{n}') * Jx;
            end
        end
    end

    function y = M_blockJacobi(~,q)
    % Solve M*y = q, where M is a block-diagonal approximation for JHJ.
        y = nan(size(q));
        offset = cache.offset;

        idx = offset(1):offset(2)-1;
        y(idx) = cache.JHJinv{1}*q(idx);
            
        idx = offset(2):offset(3)-1;
        y(idx) = tens2vec(cache.JHJinv{2}*reshape(q(idx),I,Q));
        
    end
   
end