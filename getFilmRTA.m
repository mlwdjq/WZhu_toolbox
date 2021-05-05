%% this function is used to calculate the reflectance and transmittance of a single layer film
% n0, n1, n2 are respective to the refractive index of the medium, film,
% and substrate

function [R, T, A, Rp, Tp, Rs, Ts] = ...
    getFilmRTA(n0,n1,n2,thickness_nm,lambda_nm,theta_rad)

th0 = theta_rad;
th1 = asin(real(n0).*sin(th0)./real(n1));
th2 = asin(real(n1).*sin(th1)./real(n2));
Rph = 2*pi./lambda_nm.*(n1.*thickness_nm./cos(th1) - n0.*thickness_nm.*tan(th1).*sin(th0));
Tph = 2*pi./lambda_nm.*(n1.*thickness_nm./cos(th1) - n2.*thickness_nm.*(tan(th0) - tan(th1)).*sin(th2));

% fresnel equation
[rs01,rp01,ts01,tp01] = fresnel(n0,n1,th0);
[rs12,rp12,ts12,tp12] = fresnel(n1,n2,th1);
ts10 = 1 - rs01;
tp10 = 1 + rp01;

% get film interference
[rs,ts] = fileInterference(rs01,ts01,ts10,rs12,ts12,Rph,Tph);
[rp,tp] = fileInterference(rp01,tp01,tp10,rp12,tp12,Rph,Tph);
Rs = abs(rs).^2;
Rp = abs(rp).^2;
Ts = abs(ts).^2;
Tp = abs(tp).^2;
R = (Rs + Rp)/2;
T = (Ts + Tp)/2;
A = 1 - R - T;

% fresnel equation
function [rs,rp,ts,tp] = fresnel(n0,n1,th0)
th1 = asin(real(n0).*sin(th0)./real(n1));
rs = (n0.*cos(th0) - n1.*cos(th1))./(n0.*cos(th0) + n1.*cos(th1));
rp = (n1.*cos(th0) - n0.*cos(th1))./(n1.*cos(th0) + n0.*cos(th1));
ts = 1 + rs;
tp = 1 - rp;

% get film interference
function [r,t]= fileInterference(r01,t01,t10,r12,t12,Rph,Tph)
r = r01 + t01.*t10.*r12.*exp(-1i*2*Rph)./(1+r01.*r12.*exp(-1i*2*Rph));
t = t01.*t12.*exp(-1i*Tph)./(1+r01.*r12.*exp(-1i*2*Tph));



