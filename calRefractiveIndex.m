% based on sellmeier equation to calculate refractive index
K1 = 0.34281944;
L1 = 5.74715892e-3;
K2 = 3.09568455e-1;
L2 = 1.92045655e-2;
K3 = 1.01840606;
L3 = 1.08760808e2;

lambda =0.587;

n = sqrt(K1*lambda^2/(lambda^2-L1) + ...
    K2*lambda^2/(lambda^2-L2) + K3*lambda^2/(lambda^2-L3) +1)