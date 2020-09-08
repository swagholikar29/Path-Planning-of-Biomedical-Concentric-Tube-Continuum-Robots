function [T] = arckinematics(c)
%ARCKINEMATICS Accepts a vector of k, s, theta and returns the corresponding
%homogenous transformation matrix T and (optionally) the arc link as a
%sequence of points

k   = c(1);
s   = c(2);
theta = c(3);


xRot = [0 -1 0 0;
        1 0 0 0;
        0 0 0 0;
        0 0 0 0];

xInp = @(k) [0 0 k 0;
             0 0 0 0;
            -k 0 0 1;
             0 0 0 0];

T = expm(xRot * theta) * expm(xInp(k) * s);
end