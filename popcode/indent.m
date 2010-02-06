function out = indent(in)
% out = INDENT(in) indent a string

tab = '    ';
nl = sprintf('\n');
out = [tab strrep(strtrim(in), nl, [nl tab])];