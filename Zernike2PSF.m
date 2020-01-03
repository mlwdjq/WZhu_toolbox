%% this function generate PSF based on Zernike coefficients
% please add ryan_toolbox
function PSF=Zernike2PSF(Zern,lambda, NA, xmin, xmax, res)
domain = pinhole(512);
wavefront = abdom(domain, Zern);% generate wavefront based on Zernike coefficients
H=exp(-1i*2*pi*wavefront);
PSF = abs(computePSF(H, lambda, NA, xmin, xmax, res)).^2;
imagesc(linspace(xmin,xmax,res),linspace(xmin,xmax,res),PSF);
