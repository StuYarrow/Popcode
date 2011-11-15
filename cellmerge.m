function out = cellmerge(varargin)
	s = size(varargin{1});
    for arg = 2 : length(varargin)
        if length(size(varargin{arg})) > length(s)
            s(end+1 : length(size(varargin{arg}))) = 1;
        end
        
        s = max([size(varargin{arg}) ; s], [], 1);
    end
	
	cellarr = cell(size(varargin));
	for arg = 1 : length(varargin)
		cellarr{arg} = padcell(varargin{arg}, s);
	end
	
    args = {@innermerge, cellarr{:}, 'UniformOutput', false};
    out = cellfun(args{:});
end


function merged = innermerge(varargin)        
    notEmpty = ~cellfun(@isempty, varargin);
    count = sum(notEmpty + 0);
    
    switch count
    case 0
        merged = [];
    case 1
        merged = varargin{notEmpty};
    otherwise
        error('cellmerge: conflicting data')
    end
end


function out = padcell(in, newSize)
    out = cell(newSize);
        
    for el = 1 : numel(in)
        [m n o p q r] = ind2sub(size(in), el);
        out(m,n,o,p,q,r) = in(m,n,o,p,q,r);
    end
end