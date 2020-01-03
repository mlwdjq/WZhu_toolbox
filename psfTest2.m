nT = 150;
nS = 2560;
dKLAZrn = 10*getKLAWavefront(1);
[x,y] = meshgrid(linspace(-1,1,nT));
[th,r] = cart2pol(x,y);
input_pha = zeros(nT);
for j=4:9 % using KLAWavefront as the incident aberrations
    afn = zgen([], j-1, 'fnr');
    input_pha = input_pha +dKLAZrn(j)*2*pi*afn(r,th);
end
input_pha = input_pha.*pinhole(nT);
E =pinhole(nT).*exp(1i*input_pha);
E =pad2(E,nS);
Ef = fftshift(fft2(E));
Ec = crop2(Ef,nT,nT);
%%
figure,imagesc(abs(Ec)), axis equal; axis off

