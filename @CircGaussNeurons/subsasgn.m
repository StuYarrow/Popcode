function obj = subsasgn(obj, index, val)
% SUBSASGN Define index assignment for CircGaussNeurons objects

switch index.type
case '.'
	obj = set(obj, index.subs, val);
case '()'
	error('Array indexing not supported by CircGaussNeurons objects')
case '{}'
    error('Cell array indexing not supported by CircGaussNeurons objects')
end