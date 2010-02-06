function n = minus(obj, num)
%
%
%

if isa(num, 'circdouble')
	b = num.val;
	m = num.modulo;
else
	b = num; 
end

if isa(obj, 'circdouble')
	a = obj.val;
	m = obj.modulo;
else
	a = obj;
end

n = mod(a - b, m);
n = circdouble(n, m);