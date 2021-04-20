%code to look at circ_std (angular and standard deviation) for uniform
%circular distributions for comparison

%create circular uniform distribution
thetas = [pi/80:pi/40:159*pi/80];
figure,
polarhistogram(thetas)
[angular_var,uniform_var] = circ_std(thetas,[],[],2);
p = circ_rtest(thetas)

thetas2 = [pi/8:pi/4:15*pi/8];
figure,
polarhistogram(thetas2)
[angular_var2,uniform_var2] = circ_std(thetas2,[],[],2);
p2 = circ_rtest(thetas2)

thetas3 = [0,2*pi/4,2*pi/2,6*pi/4];
figure,
polarhistogram(thetas3)
[angular_var3,uniform_var3] = circ_std(thetas3,[],[],2);
p3 = circ_rtest(thetas3)