function out = padcell(in, newSize)
    out = cell(newSize);
        
    for el = 1 : numel(in)
        [m n] = ind2sub(size(in), el);
        out(m,n) = in(m,n);
    end
end