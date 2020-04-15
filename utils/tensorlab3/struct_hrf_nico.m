function [x,state] = struct_hrf(z,task,t)
%STRUCT_HRF Haemodynamic response function.
%   [x,state] = struct_hrf(z,t,[]) computes the haemodynamical response 
%   function (HRF) x at time points in t, by modulating the model waveform 
%   with the timing parameters specified in z, as in the SPM toolbox. 
%   The structure state stores information which is reused in computing 
%   the right and left Jacobian-vector products.
%
%   struct_hrf(z,t,task) computes the right or left Jacobian-vector
%   product of this transformation, depending on the structure task. Use
%   the structure state and add the field 'r' of the same shape as z or the
%   field 'l' of the same shape as x to obtain the structure task for
%   computing the right and left Jacobian-vector products
%   
%      (dF(:)/dz(:).')*task.r(:) and
%      (dF(:)/dz(:).')'*task.l(:) + conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%   
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial
%   derivative which treats conj(z) (z) as constant. The output has the
%   same shape as x or z for the right and left Jacobian-vector products,
%   respectively.

%   REMARK: A straightforward extension is to include more gamma pdfs
%   (enlarging the dimension of the basis set), note that it might be
%   better to include these in different block terms in the decomposition
%   to also allow different scaling over voxels of the basis functions. 

%   Authors: Simon Van Eyndhoven (Simon.VanEyndhoven@esat.kuleuven.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

if nargin < 3, error('Provide time samples at which HRF should be computed'); end

if any(t<0), error('Provide only positive time samples'); end
if any(z<=0)  error('Gamma density parameters should be strictly positive'); end

if isempty(task) || (isempty(task.l) && isempty(task.r))
    % t gives the time in e.g. seconds --> convert to samples, starting at
    % zero (cfr. 'spm_hrf.m')
%     u = linspace(0,t(end)-t(1),ceil(numel(t)*(t(2)-t(1))/dt)); % upsampled time vector
    u = t(:);
    if u(1) == 0 
        u = u + (u(2)-u(1))/1000; % to avoid infinite values
    end
    % Conversion to reals needed for testing. (can be removed.)
    z = real(z);
    
    g1 = gamma_pdf(u,z(1),z(3));
    g2 = gamma_pdf(u,z(2),z(4));
%     x = g1 - g2;
    x = ((z(3)^z(1))*(u.^(z(1)-1)).*exp(-z(3)*u))/gamma(z(1)) - ((z(4)^z(2))*(u.^(z(2)-1)).*exp(-z(4)*u))/gamma(z(2));
    
    % compute the partial derivatives
    der1 =   g1.*(log(z(3)) + log(u) - psi(0,z(1)));
    der2 = - g2.*(log(z(4)) + log(u) - psi(0,z(2)));
    der3 =   g1.*(z(1)/z(3) - u);
    der4 = - g2.*(z(2)/z(4) - u);
    state.deriv = [der1(:), der2(:), der3(:), der4(:) ];
elseif ~isempty(task.r)
    x = task.deriv*(task.r.');
    state = [];
elseif ~isempty(task.l)
    x = conj(task.l'*task.deriv);
    state = [];
end

end

function gpdf = gamma_pdf(x,h,l)
    % probability density function of the gamma distribution, 
    % cfr. 'spm_Gpdf.m', assuming that h and l are scalar parameters and
    % that x does not contain negative values
    gpdf = ((l^h)*(x.^(h-1)).*exp(-l*x))/gamma(h);
end
