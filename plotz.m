function varargout = plotz(varargin)
switch length(varargin)
    case 2
        h = plot(varargin{1},varargin{2});
    case 3
        h = plot(varargin{1},varargin{2},varargin{3});
    case 4
        h = plot(varargin{1},varargin{2});
        xlabel(varargin{3});ylabel(varargin{4});
    case 5
        h = plot(varargin{1},varargin{2},varargin{3});
        xlabel(varargin{4});ylabel(varargin{5});
end
set(h,'lineWidth',2);
if ismac
    set(gca,'fontSize',16);
else
    set(gca,'fontSize',14);
end
varargout{1} = h;