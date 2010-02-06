classdef CircDouble < double

methods
    
end
function n = circdouble(varargin)
%
%
%

switch nargin
case 0
	n.val = [];
	n.modulo = [];
case 1
	if isa(varargin{1}, 'circdouble')
		n = varargin{1};
	else
		error('not a valid circdouble object')
	end
case 2
	n.val = mod(varargin{1}, varargin{2});
	i = find(n.val > (0.5 * varargin{2}));
	n.val(i) = n.val(i) - varargin{2};
	i = find(n.val < (-0.5 * varargin{2}));
	n.val(i) = n.val(i) + varargin{2};

	n.modulo = varargin{2};
otherwise
	error('wrong number of arguments')
end

n = class(n, 'circdouble');