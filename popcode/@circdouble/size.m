function s = size(varargin)
%
%
%

obj = varargin{1};
args = varargin;
args{1} = obj.val;

s = size(args{:});