function out = arraymerge(varargin)
    s = size(varargin{1});
    
    for arg = 1 : length(varargin)
        if length(size(varargin{arg})) > length(s)
            s(end+1 : length(size(varargin{arg}))) = 1;
        end
        
        s = max([size(varargin{arg}) ; s], [], 1);
    end
    
    d = length(s) + 1;
    
    arrays = {};
    for arg = 1 : length(varargin)
        sPre = size(varargin{arg});
        
        if length(s) > length(sPre)
            sPre(end+1 : length(s)) = 1;
        end
        
        arrays{arg} = padarray(varargin{arg}, s - sPre, 0, 'post');
    end
    
    stack = cat(d, arrays{:});
    out = sum(stack, d);
    check = max(abs(stack), [], d);
    
    if sum(abs(out) ~= check)
        error('arraymerge: conflicting data')
    end
end
