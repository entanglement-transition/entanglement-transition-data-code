function point = random_point_in_cylinder(R,h)
        
theta = 2 * pi * rand();
r = R*sqrt(rand());
point = [r*cos(theta), r*sin(theta), h*rand()];

end