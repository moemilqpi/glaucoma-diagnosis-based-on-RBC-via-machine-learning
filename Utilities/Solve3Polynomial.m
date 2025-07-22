function root = Solve3Polynomial(a, b, c, d)
% Refenrence: 范盛金. 一元三次方程的新求根公式与新判别法[J]. 海南师范学院学报, 1989, 2(2):91-98.

% Discriminating auxiliary variables
Delta1 = b^2 - 3*a*c;
Delta2 = b*c - 9*a*d;
Delta3 = c^2 - 3*b*d;

% Numerical accuracy adjustment
if abs(Delta1) < 1e-14, Delta1 = 0; end
if abs(Delta2) < 1e-14, Delta2 = 0; end
if abs(Delta3) < 1e-14, Delta3 = 0; end

Discriminant = Delta2^2 - 4*Delta1*Delta3;
if abs(Discriminant) < 1e-14, Discriminant = 0; end

% Situation One: Three equal real roots
if (Delta1 == 0) && (Delta2 == 0)
    x1 = -c / b;
    x2 = x1;
    x3 = x1;

% Situation Two: There exists one real root with two conjugate complex roots
elseif Discriminant > 0
    temp1 = Delta1*b + 1.5*a*(-Delta2 + sqrt(Discriminant));
    temp2 = Delta1*b + 1.5*a*(-Delta2 - sqrt(Discriminant));
    
    u = nthroot(temp1, 3);
    v = nthroot(temp2, 3);
    
    x1 = (-b - u - v) / (3*a);
    
    real_part = (-b + 0.5*(u + v)) / (3*a);
    imag_part = (sqrt(3)/2)*(u - v) / (3*a);
    
    x2 = complex(real_part, imag_part);
    x3 = complex(real_part, -imag_part);

% Situation Three: A multiple root and a single root
elseif Discriminant == 0 && (Delta1 ~= 0) && (Delta2 ~= 0)
    K = Delta2 / Delta1;
    K = round(K, 14);
    x1 = -b/a + K;
    x2 = -0.5*K;
    x3 = x2;

% Situation Four: Three different real roots
else
    sqrtDelta1 = sqrt(Delta1);
    numerator = Delta1*b - 1.5*a*Delta2;
    denominator = Delta1*sqrtDelta1;
    costheta = numerator / denominator;
    
    theta = acos(costheta);
    
    x1 = (-b - 2*sqrtDelta1*cos(theta/3)) / (3*a);
    x2 = (-b + sqrtDelta1*(cos(theta/3) + sqrt(3)*sin(theta/3))) / (3*a);
    x3 = (-b + sqrtDelta1*(cos(theta/3) - sqrt(3)*sin(theta/3))) / (3*a);
end

% Output the third root
root = x3;
end
