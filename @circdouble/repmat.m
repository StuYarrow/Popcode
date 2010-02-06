function mat = repmat(varargin)
%
%
%

obj = varargin{1};
args = varargin;
args{1} = obj.val;

mat = circdouble(repmat(args{:}), obj.modulo);