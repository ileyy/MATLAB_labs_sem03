clc
clear variables
close all

%1.3
X = randn(5);
fprintf("Исходная матрица:")
A = X * X.'

%2.2
fprintf("Встроенный метод LU разложения")
[L1, U1] = lu(A)
fprintf("Реализованный метод LU разложения")
[L2, U2] = sem_02_lu(A)
A
RES = L2 * U2

%3.2
fprintf("Встроенный метод LL разложения")
L1 = chol(A, "lower")
fprintf("Реализованный метод LL разложения")
L2 = sem_02_ll(A)
A
RES = L2 * L2.'

%4.2
fprintf("Встроенный метод QR разложения")
[Q, R] = qr(A)
fprintf("Реализованный метод QR разложения")
name = 'hausholder'
[Q, R] = sem_02_qr(A, name)
A
RES = Q * R

%2.1
function [L, U] = sem_02_lu(A)
    [n, m] = size(A);
    U = zeros(n, m);
    L = zeros(n, m);
    for i = 1 : n
        for j = i : n
            summ = 0;
            for k = 1 : i-1
                summ = summ + L(i, k) * U(k, j);
            end
            U(i, j) = A(i, j) - summ;
        end
        L(i, i) = 1;
        for j = i+1 : n
            summ = 0;
            for k = 1:(i-1)
                summ = summ + L(j, k) * U(k, i);
            end
            L(j, i) = (A(j, i) - summ) / U(i, i);
        end
    end
end

%3.1
function L = sem_02_ll(A)
    [n, n1] = size(A);
    if n ~= n1
        fprintf("ERROR. NOT SQUARE MATRIX.\n")
        return
    end

    if ~(all(eig(A)) == ones(n))
        fprintf("Matrix must be positive definite.\n")
        return
    end

    L = zeros(n);
    for i = 1 : n
        for j = 1 : i - 1
            summ = 0;
            for k = 1 : j - 1
                summ = summ + L(i, k) * L(j, k).';
            end
            L(i, j) = (1 / L(j , j)) * (A(i, j) - summ);
        end
        summ = 0;
        for k = 1 : i - 1
            summ = summ + L(i, k) * L(i, k).';
        end
        L(i, i) = sqrt(A(i, i) - summ);
    end
end

%4.1
function [Q, R] = sem_02_qr(A, name)
    [m, n] = size(A);
    if m ~= n
        error('Матрица должна быть квадратной');
    end
    
    switch name
        case 'gramm-schmidt'
            [Q, R] = gramm_schmidt(A);
        
        case 'gramm-schmidt modified'
            [Q, R] = gramm_schmidt_modified(A);
        
        case 'hausholder'
            [Q, R] = hausholder(A);
        
        case 'givens'
            [Q, R] = givens(A);
        
        otherwise
            error('Неизвестный метод: %s', name);
    end
end

function [Q, R] = gramm_schmidt(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    for j = 1 : n
        v = A(:, j);
        for i = 1 : j - 1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end

function [Q, R] = gramm_schmidt_modified(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    V = A;
    for i = 1 : n
        R(i, i) = norm(V(:, i));
        Q(:, i) = V(:, i) / R(i, i);
        for j = i + 1 : n
            R(i, j) = Q(:, i)' * V(:, j);
            V(:, j) = V(:, j) - R(i, j) * Q(:, i);
        end
    end
end

function [Q, R] = hausholder(A)
    [m, n] = size(A);
    Q = eye(m);
    R = A;
    for k = 1 : n
        x = R(k : m, k);
        e1 = zeros(length(x), 1); 
        e1(1) = 1;
        v = sign(x(1)) * norm(x) * e1 + x;
        v = v / norm(v);
        R(k : m, k : n) = R(k : m, k : n) - 2 * (v * (v' * R(k : m, k : n)));
        Q(k : m, :) = Q(k : m, :) - 2 * (v * (v' * Q(k : m, :)));
    end
    Q = Q';
end

function [Q, R] = givens(A)
    [m, n] = size(A);
    Q = eye(m);
    R = A;
    for i = 1 : n
        for j = i + 1 : m
            if R(j, i) ~= 0
                r = sqrt(R(i, i)^2 + R(j, i)^2);
                c = R(i, i) / r;
                s = -R(j, i) / r;
                G = eye(m);
                G([i j], [i j]) = [c -s; s c];
                R = G' * R;
                Q = Q * G;
            end
        end
    end
end
