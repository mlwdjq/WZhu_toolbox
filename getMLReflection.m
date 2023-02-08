function [Rs, Rp, rs_12, rp_12, ts_01, tp_01] = getMLReflection(taML,lambda_nm,AOI_deg)

n_vacuum = 1;
AOI_rad = (AOI_deg(:))'/180*pi;
nLayer = size(taML,1);
nAOI = size(AOI_rad,2);
AOIs = zeros(nLayer+1,nAOI);
AOIs(1,:) = AOI_rad;
AOIs(2,:) = asin(sin(AOIs(1,:)).*(n_vacuum)./(taML.RefractiveIndex(1)));
for j = 3:nLayer
    AOIs(j,:) = asin(sin(AOIs(j-1,:)).*(taML.RefractiveIndex(j-2))./(taML.RefractiveIndex(j-1)));
end

if nLayer > 1
   AOIs(nLayer+1,:) = asin(sin(AOIs(nLayer,:)).*taML.RefractiveIndex(nLayer-1)./(taML.RefractiveIndex(nLayer)));
end
rs_01 = zeros(nLayer,nAOI);
rp_01 = zeros(nLayer,nAOI);
rs_12 = zeros(nLayer,nAOI);
rp_12 = zeros(nLayer,nAOI);
ts_01 = zeros(nLayer,nAOI);
tp_01 = zeros(nLayer,nAOI);

%% reflectivity calculation based on Parratt's exact recursive method
% bottom layer reflection
[rs_12(1,:), rp_12(1,:)] = fresnel(taML.RefractiveIndex(nLayer),n_vacuum,AOIs(nLayer+1,:));
for j = 1:nLayer-1
    [rs_01(j,:),rp_01(j,:),ts_01(j,:),tp_01(j,:)] = fresnel(taML.RefractiveIndex(nLayer-j),taML.RefractiveIndex(nLayer+1-j),AOIs(nLayer+1-j,:));
    %total reflection from the bottom layer
    rs_12(j+1,:) = slab(rs_01(j,:),rs_12(j,:),taML.Thickness(nLayer+1-j),...
        lambda_nm,taML.RefractiveIndex(nLayer-j),taML.RefractiveIndex(nLayer+1-j),AOIs(nLayer+1-j));
    rp_12(j+1,:) = slab(rp_01(j,:),rp_12(j,:),taML.Thickness(nLayer+1-j),...
        lambda_nm,taML.RefractiveIndex(nLayer-j),taML.RefractiveIndex(nLayer+1-j),AOIs(nLayer+1-j));
end

[rs_01(nLayer,:),rp_01(nLayer,:),ts_01(nLayer,:),tp_01(nLayer,:)] = fresnel(n_vacuum,taML.RefractiveIndex(1),ALOs(1,:));

% total reflection from the cap layer
Rs = slab(rs_01(nLayer,:),rs_12(nLayer,:),taML.Thickness(1),lambda_nm,n_vacuum,taML.RefractiveIndex(1),AOIs(1,:));
Rp = slab(rp_01(nLayer,:),rp_12(nLayer,:),taML.Thickness(1),lambda_nm,n_vacuum,taML.RefractiveIndex(1),AOIs(1,:));





