function varargout = subsref(obj, index)
%
%
%

switch index.type
case '.'
    switch index.subs
    case 'val'
    	varargout = {obj.val};
    case 'modulo'
        varargout = {obj.modulo};
    otherwise
        error('Invalid field name')
    end

case '()'
	varargout = {obj.val(index.subs{:})};
	%error('Array indexing not supported by circdouble objects')
case '{}'
    error('Cell array indexing not supported by circdouble objects')
end