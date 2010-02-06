function obj = widthadapt(obj, width, amnt, centre)
%
%
%

obj = set(obj, 'width', obj.width .* (1 - amnt .* exp(-(1.0 ./ deg2rad(width).^2) .* (1 - cosd(double(obj.Neurons.preferredStimulus - centre))))));