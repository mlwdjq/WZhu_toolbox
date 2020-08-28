function [E,n,m] = elliptical_zernike5(r,phi,N, b, varargin)% ordering, orthonormal, orientation)
%%
% This function generates polynomials orthogonal (orthonormal) across an
% elliptical pupil. These are equivalent to Zernike polynomials on a
% circular pupil. The method is based on Refs [1] and [2].
% [1] JOSA A 24,7 (2007)
% [2] Optics Letters 32, 1 (2007)
% Note: eq. (14) in Ref. [1] has a typo. It should be normalized by the
% area of the pupil, as in Ref. [2].

% Inputs:
% r = radius-vector (2D matrix)
% phi = azimuthal angle wrt the x-axis (2D matrix)
% Note: The best way to generate r and phi is from a uniform cartesian grid
% [x,y] = meshgrid(-1:dx:1);
% [phi, r] = cart2pol(x,y);
% N = polynomial order; if numel(N) == 2, ordering is ignored.
% If numel(N) == 1, ordering is required
% b = minor(y)/major(x) axis. 0<b<=1
% ordering = index ordering type for Zernike polynomials (prolith or noll)
% orthonormal (optional) = 1 to normalize polynomials,
%                          0 (default) to leave unnormalized
% orientation (optional) = 'h' puts the major axis along x,
%                          'v' (default) puts the major axis along y,

% Outputs:
% E = polynomial values on a grid defined by r and phi
% n = radial polynomial order
% msgn = azimuthal polynomial order


% Author: Dmitriy Zusin, Cube 5.10H-167. December 26, 2018.

%% Parse the inputs
p = inputParser;
addRequired(p, 'r');
addRequired(p, 'phi');
addRequired(p, 'N');
addRequired(p, 'b');
addOptional(p,'ordering', '', @ischar);
addOptional(p, 'orthonormal', 0, @isnumeric);
addOptional(p, 'orientation', 'v', @ischar);
parse(p, r,phi,N, b, varargin{:});

ordering = p.Results.ordering;
orthonormal = p.Results.orthonormal;
orientation = p.Results.orientation;

%% Index conversion
% Convert to Noll's index convention, which is used internally, according
% to Ref. [1]
Nn = zeros(size(N,1), 1);
n = Nn;
m = Nn;

for i = 1:numel(Nn)
    % Index conversion: (extract the radial (n) and azimuthal (m) polynomial
    % orders)
    [n(i), m(i)] = ind2nm(N(i,:), ordering);
    Nn(i) = ind2noll(n(i), m(i));
end

N = max(Nn);

% query cartesian points
[xq,yq] = pol2cart(phi, r);
x = xq; y = yq;

%% Grid generation
% create a fine cartesian grid. A fine grid is necessary for an accurate
% computation of the polynomial conversion matrix elements, which involves
% integration over the pupil area. A coarse grid results in a large
% numerical integration error.

% Create a uniform square grid from -1 to 1 if [x,y] is not uniform
dx = diff(x,1,2);
dx = dx(dx>0);
dx = min(dx(:));
dy = diff(y);
dy = dy(dy>0);
dy = min(dy(:));

% 0.02 is the largest grid step size allowed. Otherwise the numerical
% integration becomes inaccurate.
Ngrid = round(2./ min([dx, dy, 0.02])) + 1;
[x,y] = meshgrid(linspace(-1,1,Ngrid));

% Prepare a [phi, r] grid for zernike5.m
[phi, r] = cart2pol(x,y);

%% Elliptical pupil
switch orientation
    case 'v'
        pupil = x.^2./b.^2 + y.^2 <= 1;
    case 'h'
        pupil = x.^2 + y.^2./b.^2 <= 1;
end

%% Circular Zernike polynomials
% Compute Zernike polynomials up to order N and store them in Z
% [1] defines orthonormal Zernike polynomials (see eqs. (1a-c)). Hence the
% sqrt(n+1) factor below.
Z = zeros(size(pupil,1), size(pupil,2), N); % allocate Z
for j = 1:N
    [Z(:,:,j), n1, m1] = zernike5(r, phi, j, 'noll');
    Z(:,:,j) = Z(:,:,j) .*sqrt(n1+1);
    if m1~=0, Z(:,:,j) = sqrt(2) .* Z(:,:,j); end
end

%% Matrix of inner products
% Allocate a matrix of inner products of Zernike polynomials
C = zeros(N, N);

for i = 1:N
    for j = i:N
        tmp = Z(:,:,i).*Z(:,:,j) .* pupil;
        C(i,j) = sum(tmp(:));
    end
end

% Note: there's a typo in eq. (14) in Ref. [1]. It should normalize by the
% pupil area, which is C(1,1).
C = C ./ C(1,1);

%% Elliptical polynomials
E = zeros(size(xq,1), size(xq,2), length(Nn));

for i = 1:length(Nn)
    
    % Conversion matrix (solve for Q from Q'*Q = C using Cholesky
    % factorization, according to Ref. [1] eqs. (12) and (13)
    M = inv(chol(C(1:Nn(i), 1:Nn(i)))');
    
    % Find the elliptical polynomial, according to eq. (11) in Ref. [1]. E is
    % orthonormal by the procedure.
    E_tmp = zeros(size(x));
    for j = 1:Nn(i)
        E_tmp = E_tmp + M(Nn(i), j) .* Z(:,:,j);
    end
    
    % Compute values of E at the query grid points
    E(:,:,i) = interp2(x,y, E_tmp, xq, yq);
    
    % Unnormalize if necessary
    if orthonormal == 0
        E(:,:,i) = E(:,:,i) ./ sqrt((n(i)+1) .* (2.*(m(i)~=0) + 1.*(m(i)==0) ));
    end
    
end

E = squeeze(E);

end

%%
% An auxiliary nested function to compute Noll's indices
function noll_ind = ind2noll(n,m)

rem = mod(n,4);
noll_ind = n*(n+1)/2 + abs(m);
if m>=0 && (rem == 2 || rem == 3)
    noll_ind = noll_ind + 1;
end

if m<=0 && (rem == 0 || rem == 1)
    noll_ind = noll_ind + 1;
end

end

%%
% An auxiliary nested function to convert a polynomial index to the radial
% and azimuthal polynomial numbers
function [n, m] = ind2nm(N, ordering)

switch numel(N)
    case 2
        n = N(1);
        m = N(2);
    case 1
        switch lower(ordering)
            case 'noll'
                n = 1:N;
                n = find( (n+1).*n/2 < N , 1, 'last' ) ;
                if isempty(n)
                    n=0;
                end
                l = n*(n+1)/2+1;
                m = ( (N-l)  + rem( N-l+n, 2) ) * (-1)^N;
            case 'prolith'
                q = floor(sqrt(N));
                if q^2 == N
                    n = (q-1)*2;
                    m = 0;
                else
                    % m + n = 2 * q
                    p = (q+1)^2 - N;
                    m = ceil( p/2 ) * (-1)^rem(p,2);
                    n = 2 * q - abs(m);
                end
            otherwise
                error('Unknown indexing scheme');
        end
    otherwise
        error('index N in zernike5(r,phi,N,ordering)is not understood')
end

end