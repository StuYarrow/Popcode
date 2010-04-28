function out = cellmerge(varargin)
	s = size(varargin{1});
    for arg = 2 : length(varargin)
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
