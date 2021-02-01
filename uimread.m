function [I,map]=uimread()
%% read image
I=0;map=0;
global old_pr
if old_pr==0
    old_pr=[];
end
[fn,pn]=uigetfile({'*.bmp;*.jpg; *.tiff;*.tif; *.gif; *.png','Image';'*.*','all'},'Load image',old_pr);
old_pr=pn;
if isequal(fn,0)||isequal(pn,0)
    return;
end
[I,map]=imread([pn fn]);