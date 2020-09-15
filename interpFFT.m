function Es = interpFFT(E,M,N)
[ysize,xsize] = size(input);
spectrum = fftshift(fft2(fftshift(E)))/sqrt(ysize*xsize);
spectrum = pad2(spectrum,M,N);
Es = ifftshift(ifft2(ifftshift(spectrum)))*sqrt(M*N);