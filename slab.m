% single layer reflection equation
function r = slab(r01,r12,d_nm,lambda_nm,n0,n1,th0)
delta = 4*pi*d_nm./lambda_nm.*sqrt(n1^2-n0^2.*sin(th0));
r = (r01+r12.*exp(-1i*delta))./(1+r01.*r12.*exp(1i*delta));
