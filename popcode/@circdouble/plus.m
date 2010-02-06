function n = plus(obj, num)
%
%
%

if isa(num, 'circdouble')
	num = num.val;
end

n = circdouble(mod(obj.val + num, obj.modulo), obj.modulo);