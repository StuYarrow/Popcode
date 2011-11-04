function h = singleEBlin(x, y, err, width, label, col, font)

le = errorbar(x, y, err);
dx = diff(get(gca,'XLim'));
errorbar_tick(le, dx/width)
set(le, 'linestyle', '-', 'color', col, 'marker', 'none')

line([x - width*0.25, x, x + width*0.25], [y y y], 'color', col)

t = text(x + width * 0.75, y, label);
set(t, 'fontsize', font, 'verticalalignment', 'middle', 'horizontalalignment', 'left', 'verticalalignment', 'middle')

h = [le t];