function out = cellmerge(varargin)
    args = {@innermerge, varargin{:}, 'UniformOutput', false}
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
