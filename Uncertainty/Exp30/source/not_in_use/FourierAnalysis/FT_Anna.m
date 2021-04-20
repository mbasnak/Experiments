function [phase,amplitude] = FT_Anna(signal)
% very drafty, assumes data to be between [-Pi,Pi}


% for testing
% signal = [12 + 0.3 * cos(x'+0.2) + 4 * cos(2*x' - 0.4 ) , ...
%            25 * cos(x') + 9 * cos(2*x'-0.9) ];

N = size(signal,1);
dx = 2 * pi/N;
x = linspace(-pi,pi-dx,N);

k = 0:(N-1);
ft_cos = 1/pi * cos(k' * x) * signal * dx;
ft_sin = 1/pi * sin(k' * x) * signal * dx;

phase = atan2(ft_sin,ft_cos);
amplitude = sqrt(ft_cos.^2 + ft_sin.^2);

end

