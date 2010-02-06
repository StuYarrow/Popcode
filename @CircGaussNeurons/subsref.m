function varargout = subsref(obj, index)
%SUBSREF Define field name indexing for Neurons objects


switch index.type
case '.'
	varargout = {get(obj, index.subs)};
case '()'
	error('Array indexing not supported by Neurons objects')
case '{}'
    error('Cell array indexing not supported by Neurons objects')
end