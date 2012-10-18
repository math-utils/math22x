function out = gauss_blur(x,y,x1,fwhm)
% Jonathan Ong, jongleberry@gmail.com. 2010
% Uses Gaussian functions to resample/smooth/convolute/blur a given
% function.
% x are the locations of the samples. It does not have to be uniform.
% y are the values of the sample.
% x1 are the output locations you want to resample to. 
% fwhm is the full width at half maximum: 
% http://en.wikipedia.org/wiki/Full_width_at_half_maximum
% fwhm can be either the same size as x1 or a scalar
% This code assumes x,y,x1 are all column vectors (one column)
% To use sparse matrix representation, all inputs must be in double
% precision
if ~strcmp(class(x),'double')
    x = double(x);
end
if ~strcmp(class(y),'double')
    y = double(y);
end
if ~strcmp(class(x1),'double')
    x1 = double(x1);
end
[x,i] = sort(x);
y = y(i);
if length(x) ~= length(y)
    error('Inputs x and y must be the same length')
end
if (length(x1) ~= length(fwhm)) && (length(fwhm) ~= 1)
    error('fwhm must be scalar or the same size as x1')
end
n = length(x);
nout = length(x1);
if length(fwhm) == 1
    fwhm = fwhm*ones(nout,1);
end
A = exp(-4*log(2)*((repmat(x,1,nout)-repmat(x1',n,1))./repmat(fwhm',n,1)).^2);
A(abs(A) < 10*eps) = 0;
A = sparse(A);
curv = zeros(size(x));
curv(2:end-1) = (x(3:end)-x(1:end-2))/2;
curv(1) = x(2)-x(1);
curv(end) = x(end)-x(end-1);
A = A.*repmat(curv,1,nout);
out = sum(sparse(diag(y))*A,1)./sum(A,1);