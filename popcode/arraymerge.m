function out = arraymerge(varargin)
    s = size(varargin{1});
    d = length(s) + 1;
    
    for arg = 1 : length(varargin)
        s = max([size(varargin{arg}) ; s], [], 1);
    end
    
    arrays = {};
    for arg = 1 : length(varargin)
        arrays{arg} = padarray(varargin{arg}, s - size(varargin{arg}), 0, 'post');
    end
    
    stack = cat(d, arrays{:});
    out = sum(stack, d);
    check = max(abs(stack), [], d);
    
    if sum(abs(out) ~= check)
        error('arraymerge: conflicting data')
    end
end
