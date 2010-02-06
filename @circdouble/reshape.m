function mat = reshape(varargin)
%
%
%

obj = varargin{1};
args = varargin;
args{1} = obj.val;
mat = circdouble(reshape(args{:}), obj.modulo);