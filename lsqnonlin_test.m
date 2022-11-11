% Testing if "trick" works for lsqnonlin
% Original Cost Function is (x - 3)^2 + (y - 5)^2
% Linear inequality: x >= y
%
% -x + y <= 0
% --> delta = -x + y, with upper bound of 0

X0 = [100,200];
X = lsqnonlin(@cost,X0,[-inf,-inf],[inf,0]);


function F = cost(vars)
    x = vars(1); delta = vars(2);
    y = x + delta;
    F = [x - 3;
         y - 5];
end