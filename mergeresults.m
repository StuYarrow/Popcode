function varargout = mergeresults(files, vars)
    
    out = cell(1,length(vars));
    
    for v = 1 : length(vars)
        varname = vars{v};
        in = cell(1,length(files));
        
        for f = 1 : length(files)
            filename = files{f};
            inCell = struct2cell(load(filename, varname));
            in(f) = inCell(1);
        end
                
        if iscell(in{1})
            out(v) = cellmerge(in{:});
        elseif isnumeric(in{1})
            out{v} = arraymerge(in{:});
        else
            error('mergeresults: not cell or numeric')
        end     
    end
    
    if nargout == 1
        varargout = {out};
    elseif nargout == length(out)
        varargout = out;
    else
        error('mergeresults: wrong number of outputs')
    end
end