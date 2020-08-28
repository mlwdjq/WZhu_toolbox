function Ecoeff = WF2EZernike(WF, b, Ne)
%%
% This function computes wavefront maps, i.e., expansion coefficients onto
% a set of elliptical polynomials.
%
% Inputs:
% WF = a 3D array of wavefronts. y and x are along dimensions 1 and 2.
%      The 2D wavefront images are stacked along the third dimension. The
%      wavefront pupil has to fill the entire grid along at least dimension
%      1 or 2. Preferably, the borders of the window should be tight around
%      the pupil. The import_intfile.m function already takes care of it.
% b = pupil ellipticity (scalar).
% Ne = highest polynomial order up to which the expansion is computed. The
%      expansion coefficients are calculated from order to Ne
%
% Outputs:
% Ecoeff = 2D matrix of expansion coefficients in mWv. The rows
%   correspond to different wavefronts in the WF stack. The number of rows
%   is, therefore, equal to size(WF,3). The columns correspond to different
%   polynomial orders that go from 2 to Ne. The number of columns is equal
%   to Ne-1.
%
% Author: Dmitriy Zusin, January 29, 2019. Cube 10H-167, Building 5.

%%

% Determine the pupil based on the wavefront data. The pupil represents the
% region with valid data, which is not necessarily strictly elliptical.
pupil_wf = ~isnan(mean(WF,3));

% pupil centroid
xc = mean(find(sum(pupil_wf) == max(sum(pupil_wf))));
yc = mean(find(sum(pupil_wf,2) == max(sum(pupil_wf,2))));

j=(2:Ne)'; % polynomial orders

% cartesian grid
N = max(size(WF(:,:,1)));
[x, y] = meshgrid(linspace(-1,1,N));

% Generate a strictly elliptical pupil on a cartesian grid. The pupil
% orientation is determined by the aspect ratio of the pupil sides.
if size(WF,1)>=size(WF,2)
    orientation = 'v';
    pupil = x.^2./b.^2 + y.^2 <= 1;
else
    orientation = 'h';
    pupil = y.^2./b.^2 + x.^2 <= 1;
end

% Centroid of the strictly elliptical pupil
xc0 = mean(find(sum(pupil) == max(sum(pupil))));
yc0 = mean(find(sum(pupil,2) == max(sum(pupil,2))));

% pad the wavefront in order to make its window square and put the center
% of the pupil at the center of the grid
if size(WF,1)>size(WF,2)
    WF = padarray(WF, [0 (xc0-xc + 1) 0], NaN, 'pre');
    WF = padarray(WF, [0 (size(x,2)-size(WF,2)) 0], NaN, 'post');
elseif size(WF,1)<size(WF,2)
    WF = padarray(WF, [(yc0-yc + 1) 0 0], NaN, 'pre');
    WF = padarray(WF, [(size(y,1)-size(WF,1)) 0 0], NaN, 'post');
end

% Compute elliptical polynomials
[phi, r] = cart2pol(x,y);
E = elliptical_zernike5(r,phi,j, b, 'prolith', 0, orientation);

% the smallest pupil is the product of pupil and pupil_wf, in case the two
% do not perfectly agree.
pupil_wf = ~isnan(mean(WF,3));
pupil = pupil .* pupil_wf;
pupil(~pupil) = NaN;

% determine the indicies of NaN's for future removal
pupil = pupil(:);
ind = find(isnan(pupil));

% reshape and normalize the elliptical polynomials. Get rid of the elements
% corresponding to NaN's in the pupil
E = reshape(E, size(E,1).*size(E,2), size(E,3));
E(ind,:) = [];
E = E ./ std(E); % this step is important. It determines the normalization

% Allocate a matrix of coefficients
Ecoeff = zeros(size(WF,3), length(j));

% Project the wavefronts onto the elliptical polynomials
for i = 1:size(WF,3)
    
    WF_i = squeeze(WF(:,:,i));
    WF_i = WF_i(:);
    WF_i(ind) = []; % remove the indices corresponding to NaN's in the pupil
    
    Ecoeff(i,:) = sum(WF_i .* E) ./ sum( E.* E ); % expansion coefficients
end