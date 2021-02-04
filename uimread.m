function [I,map]=uimread(multiSelection)
% load images
if nargin == 0
    multiSelection = 'off';
end
if ~strcmp(multiSelection,'on')
    multiSelection = 'off';
end
I=0;map=0;
global old_pr
if old_pr==0
    old_pr=[];
end
[fn,pn]=uigetfile({'*.bmp;*.jpg; *.tiff;*.tif; *.gif; *.png','Image';'*.*','all'},'Load image',old_pr,'MultiSelect',multiSelection);
old_pr=pn;
if isequal(fn,0)||isequal(pn,0)
    return;
end
if iscell(fn)
    [~,N]=size(fn);
    [Is,map]=imread([pn fn{1}]);
    I = repmat(Is(:,:,1),1,1,N);
    for i=2:N
        [Is,map]=imread([pn fn{i}]);
        I(:,:,i)=Is(:,:,1);
    end
else
    [Is,map]=imread([pn fn]);
    I=Is(:,:,1);
end
