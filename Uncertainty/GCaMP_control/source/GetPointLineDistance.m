function distance = GetPointLineDistance(x3,y3,x1,y1,x2,y2)
% Get the distance from a point (x3, y3) to
% a line defined by two points (x1, y1) and (x2, y2);

	% Find the numerator for our point-to-line distance formula.
	numerator = abs((x2 - x1) * (y1 - y3) - (x1 - x3) * (y2 - y1));
	
	% Find the denominator for our point-to-line distance formula.
	denominator = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
	
	% Compute the distance.
	distance = numerator ./ denominator;

end