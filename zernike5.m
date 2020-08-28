function  [Z,n,msgn] = zernike5(r,phi,N,ordering)
% [Z,n,msgn] = zernike5(r,phi,N,ordering)
% x,y are cartesian coordinates of points in the unit circle
% x and y must have the same size. 
% On output, size(Z,1) = numel(x) = numel(y); 
% Z(:,n) is the value of the n-th Zernike polynomial at [x(:),y(:)]
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
        end
    otherwise
        error('index N in zernike5(r,phi,N,ordering)is not understood')
end

if m < 0
    Z = -sin( m * phi ); 
elseif m > 0
    Z = cos( m * phi ); 
else
    Z = ones(size(r)); 
end
%fprintf('%d %d\n',n,m)
if isempty(r)
    return
end
msgn= m; 
m   = abs(m); 
npm = ( n + m )/2; 
nmm = ( n - m )/2; 

assert( round(npm) == npm ); 
assert( round(nmm) == nmm ); 

fac = prod( (npm+1) : n ) / factorial( nmm ); 
Rnm = 0; 
for k=0:nmm
   %fprintf('%g ',fac)
   Rnm = Rnm +  fac * r.^(n-2*k); 
   fac = - fac * ( npm - k ) * ( nmm - k ) / ( (k+1) * ( n - k ) );
end
Z = Z .* Rnm; 
%fprintf('\n')
end