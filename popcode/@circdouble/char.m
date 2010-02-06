function dispStr = char(obj)
%
%
%

dispStr =  strvcat(['(' num2str(obj.modulo) ')'], num2str(obj.val));