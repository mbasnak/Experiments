function [x_out, y_out] = removeWrappedLines(x, y)

% inserts a NaN between points where a wrapping occurs, so that when
% plotted the wrap doesn't show up

% takes in degrees, looks for jumps of 180 degrees
y_out = [];
x_out = [];

ind_jump = find(abs(diff(y))>100);

if ~isempty(ind_jump)
    for i = 1:length(ind_jump)
        index = ind_jump(i);
        if i == 1
            index_before = 1;
        else
            index_before = ind_jump(i-1)+1;
        end
        x_out = [x_out x(index_before:index)' NaN];
        y_out = [y_out y(index_before:index)' NaN];
    end

    index_before = ind_jump(end)+1;
    if index_before < length(x)
        x_out = [x_out x(index_before:end)'];
        y_out = [y_out y(index_before:end)'];
    end
else
    x_out = x;
    y_out = y;
end

end