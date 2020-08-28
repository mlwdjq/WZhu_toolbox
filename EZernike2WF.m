function WF = EZernike2WF( Ecoeff, b, N, varargin)
%Combine Zernike coefficients back into a 2D wavefront map
%Input: Zc -- Zernike coefficient vector in rms convention
%       sz -- size of the wavefront map
%Output: WF -- 2D map in waves

%%
% This function combines elliptical Zernike coefficients back into a 2D
% wavefront map. It can do so at once for mutliple sets of expansion
% coefficients
%
% b = pupil ellipticity (scalar).
%
% Inputs:
% Ecoeff = 2D matrix of expansion coefficients in mWv. The rows
%   correspond to different wavefronts in the WF stack. The number of rows
%   is, therefore, equal to size(WF,3). The columns correspond to different
%   polynomial orders that go from 2 to size(Ecoeff,2)+1.
% b = ellipticity of the pupil
% N = number of grid points. If N is less than 801, it will automatically
%   default to 801.
% orientation (optional) = 'h' puts the major axis along x,
%                          'v' (default) puts the major axis along y
%
% Outputs:
% WF = a 3D array of wavefronts. y and x are along dimensions 1 and 2.
%      The 2D wavefront images are stacked along the third dimension.
%
%
% Author: Dmitriy Zusin, January 29, 2019. Cube 10H-167, Building 5.

%% Parse the inputs
p = inputParser;
addRequired(p, 'Ecoef');
addRequired(p, 'b');
addRequired(p, 'N');
addOptional(p, 'orientation', 'v', @ischar);
parse(p, Ecoeff, b, N, varargin{:});

orientation = p.Results.orientation;
%% Spatial grid
[x,y] = meshgrid(linspace(-1,1,N));
[phi, r] = cart2pol(x,y);

%% Elliptical pupil
switch orientation
    case 'v'
        pupil = x.^2./b.^2 + y.^2 <= 1;
    case 'h'
        pupil = x.^2 + y.^2./b.^2 <= 1;
end

% polynomial orders
j = 1:size(Ecoeff,2);
j = (j+1)';

% generate elliptical polynomials
E = elliptical_zernike5(r,phi,j, b, 'prolith', 0, orientation);

% compute the normalization factor
pupil = ~isnan(sum(E,3)) & pupil;
norm = reshape(E, size(E,1).*size(E,2), size(E,3));
norm = std(norm(pupil(:),:));
norm = reshape(norm, 1, 1, length(norm));

% allocate the wavefront array
WF = zeros(size(E,1), size(E,2), size(Ecoeff,1));

% Combine the polynomials into a wavefront map
for i = 1:size(Ecoeff, 1)
    Ecoeff_tmp = reshape(Ecoeff(i,:), 1, 1, size(Ecoeff, 2));
    WF_tmp = sum(Ecoeff_tmp ./ norm .* E .* pupil, 3);
    WF_tmp(~pupil) = NaN;
    WF(:,:,i) = WF_tmp;
end

% remove singleton dimensions
WF = squeeze(WF);