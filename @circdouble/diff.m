function n = diff(varargin)
%
%
%

args = varargin;
obj = varargin{1};
args{1} = obj.val;

n = circdouble(diff(args{:}), obj.modulo);