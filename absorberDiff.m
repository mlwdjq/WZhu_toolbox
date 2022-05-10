function Ej = absorberDiff(Ei,Mn,dxy,dz,lambda,theta,phi,direction)

[ny,nx,nz] = size(Mn);

[X,Y] = meshgrid((1:nx)*dxy,(1:ny)*dxy);

% transmission of each slice
T = exp(1i*2*pi/lambda*Mn*dz*cos(theta));

n_vacuum = 1;

E0 = exp(1i*2*pi/lambda*n_vacuum*sin(theta)*cos(phi)*X).*exp(1i*2*pi/lambda*n_vacuum*sin(theta)*sin(phi)*Y);

% propagation
Ej = Ei;
fx = (-ceil((nx-1)/2):floor((nx-1)/2))/(nx*dxy) + sin(theta)*cos(phi)/lambda;
fy = (-ceil((ny-1)/2):floor((ny-1)/2))/(ny*dxy) + sin(theta)*sin(phi)/lambda;
[FX,FY] = meshgrid(fx,fy);

for j = 1:nz
    if direction ==1
        Tj = T(:,:,nz+1-j);
    else
        Tj = T(:,:,j);
    end
    Ej = Ej.*Tj;
    fEj = fftshift(fft2(Ej.*conj(E0))); impose pseudo-periodicity for arbitrary incident angle
    fEj = fEj.*exp(1i*2*pi/lambda*dz*(sqrt(1-lambda^2*FX.^2-lambda^2*FY.^2)-1));
    Ej = ifft2(fftshift(fEj)).*E0;
end
