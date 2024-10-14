clc
clear variables 
close all

A = [1 -0.2589 -0.3093,
    -0.2589 1 -0.2705,
    -0.3093 -0.2705 1];

b = [1,
    1,
    1];

x = [2.2873,
    2.2162,
    2.3068];

x1 = sem_03_matrix(A, b);
x2 = sem_03_kramer(A, b);
x3 = sem_03_gauss(A, b);
x4 = sem_03_gauss_jordan(A, b);

fprintf('Size of matrix A: %d x %d\n\n', size(A));
fprintf('Size of matrix b: %d x %d\n\n', size(b));
disp(table(A, b, x, x1, x2, x3, x4));

%2.1
function [x, ok] = sem_03_matrix(A, b)
    if det(A) == 0
        ok = false;
        x = zeros(size(b));
        return;
    end
    
    try
        x = inv(A) * b;
        ok = true;
    catch
        ok = false;
        x = zeros(size(b));
    end
end

function [x, ok] = sem_03_kramer(A, b)
    n = size(A, 1);
    if det(A) == 0 || size(A, 1) ~= size(A, 2)
        ok = false;
        x = zeros(size(b));
        return;
    end
    
    x = zeros(n, 1);
    detA = det(A);
    
    try
        for i = 1 : n
            Ai = A;
            Ai(:, i) = b;
            x(i) = det(Ai) / detA;
        end
        ok = true;
    catch
        ok = false;
        x = zeros(n, 1);
    end
end

function [x, ok] = sem_03_gauss(A, b)
    n = length(b);
    x = zeros(n, 1);
    ok = true;

    for i = 1 : n
        if A(i, i) == 0
            ok = false;
            return;
        end
    end

    for i = 1 : n
        for j = i + 1 : n
            k = A(j, i) / A(i, i);
            A(j, :) = A(j, :) - (A(i, :) * k);
            b(j) = b(j) - k * b(i);
        end
    end
    for i = n : -1 : 1
        x(i) = b(i);
        for k = i + 1 : n
            x(i) = x(i) - A(i, k) * x(k);
        end
        x(i) = x(i) / A(i, i);
    end
end

function [x, ok] = sem_03_gauss_jordan(A, b)
    n = length(b);
    aug = [A, b];
    
    for i = 1 : n
        if aug(i, i) == 0
            ok = false;
            x = zeros(n, 1);
            return;
        end
        
        aug(i, :) = aug(i, :) / aug(i, i);
        
        for j = 1 : n
            if j ~= i
                factor = aug(j, i);
                aug(j, :) = aug(j, :) - factor * aug(i, :);
            end
        end
    end
    
    x = aug(:, end);
    ok = true;
end
