% this function is used to rotate/flip the zernike aberrations to rotated
% zernike basis
% default rot = 1, turn 90 deg clockwise 
% default flip = 1; no flip

function [dZern, A] = rotateZernike(dZern, rot, flip, nZern)
N = 100;
if nargin ==3
    nZern = 37;
elseif nargin ==2
    nZern = 37;
    flip = 1;
elseif nargin ==1
    nZern = 37;
    flip = 1;
    rot = 1;
end
mask= pinhole(N);
% generate rotation matrix
A = zeros(nZern);
for i = 1:nZern
    abvec = zeros(nZern,1);
    abvec(i) = 1;
    wave = abdom(mask, abvec);
    A(:,i) = zerndecomp(flip*rot90(wave,rot), nZern-1, mask);
end
dZern = A * dZern(:);