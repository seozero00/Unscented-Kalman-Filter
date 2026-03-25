function [xm, xcov] = UT(Xi, W, noiseCov)
[n,kmax] = size(Xi);

xm = 0;
for k=1:kmax
    xm = xm + W(k)*Xi(:, k); % xm = sigma(Wi*xi) (i=1:2n+1)
end

xcov = zeros(n,n); % xcov = Px = sigma(Wi*(xi-xm)*(xi-xm)') (i=1:2n+1)
for k=1:kmax
    xcov = xcov + W(k)*(Xi(:, k) - xm )*(Xi(:, k) - xm)';
end
xcov = xcov + noiseCov;