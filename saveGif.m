function saveGif(h,filename,delay)
if nargin==0
    h=gcf;
    filename='untitled.gif';
    delay=0.05;
elseif nargin==1
    filename='untitled.gif';
    delay=0.05;
elseif nargin==2
    delay=0.05;
end
frame = getframe(h);
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
if ~exist(filename)
    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delay);
else
    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay);
end


