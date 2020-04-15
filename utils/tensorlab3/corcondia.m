function [consistency] = corcondia(T,U)
%CORCONDIA Core consistency diagnostic for CPD.
%   CONSISTENCY = CORCONDIA(T,U) computes the core consistency diagnostic for a given tensor T
%   and its approximation CPD, given as a cell of factor matrices U. 
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/04/06   NV      Initial version

    N = getorder(T);
    R = size(U{1},2);
    if numel(U) > N
        % Remove additional 'one'-modes
        nrm  = prod(vertcat(U{N+1:end}));
        U{1} = bsxfun(@times, U{1}, nrm);
        U    = U(1:N);
    end 
    % Normalize factors for improved stability
    nrm = cellfun(@(u) sqrt(sum(abs(u).^2)), U, 'UniformOutput', false);
    U   = cellfun(@(u,n) bsxfun(@rdivide, u, n), U, nrm, 'UniformOutput', false);
    nrm = prod(vertcat(nrm{:}),1);
    % Compute approximating core tensor
    U = cellfun(@pinv, U, 'UniformOutput', false);
    Sest = reshape(mtkronprod(T, U, 0, 'H'), ones(1,N)*R);
    % Compute metric
    nrmdiag = sum(abs(Sest(linspace(1,R^N,R))).^2);
    Sest(linspace(1,R^N,R)) = Sest(linspace(1,R^N,R)) - nrm;
    consistency = 100*(1 - frob(Sest,'squared') / nrmdiag);
end
