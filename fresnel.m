% fresnel equation
function [rs,rp,ts,tp] = fresnel(n0,n1,th0)
th1 = asin((n0).*sin(th0)./(n1));
rs = (n0.*cos(th0) - n1.*cos(th1))./(n0.*cos(th0) + n1.*cos(th1));
rp = (n1.*cos(th0) - n0.*cos(th1))./(n1.*cos(th0) + n0.*cos(th1));
ts = 1 + rs;
tp = 1 - rp;