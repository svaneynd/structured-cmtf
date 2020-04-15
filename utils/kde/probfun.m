function pdf=probfun(x,w,mu,Sig)

% Reference:
%  Kernel density estimation via diffusion
%  Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
%  Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
%  https://nl.mathworks.com/matlabcentral/fileexchange/58312-kernel-density-estimator-for-high-dimensions?s_tid=prof_contriblnk


[gam,d]=size(mu);
pdf=0;
for k=1:gam
    L=chol(Sig(:,:,k));s=diag(L);
    logpdf=-.5*sum(( bsxfun(@minus,x,mu(k,:))/L).^2,2)+log(w(k))...
        -sum(log(s))-d*log(2*pi)/2;
    pdf=pdf+exp(logpdf);
end
end