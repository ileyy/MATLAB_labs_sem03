clc
clear variables
close all

%1.2
N_a = 8;
N_b = 10 * N_a;

x1 = linspace(0, 1, N_a);
x0 = linspace(0, 1, N_b);

%1.3
f0 = func(x0);
f1 = func(x1);

%3.1
f_nn = sem_04_nn(x1, f1, x0);
f_linear = sem_04_linear(x1, f1, x0);
f_lagrange = sem_04_lagrange(x1, f1, x0);
f_newton_forward = sem_04_newton_forward(x1, f1, x0, 3);
f_newton_backward = sem_04_newton_backward(x1, f1, x0, 3);

%3.2
figure;

subplot(2, 1, 1);
hold on;
plot(x0, f0, 'k', 'LineWidth', 1.5, 'DisplayName', 'Аналитическая функция');

stem(x1, f1, 'r', 'filled', 'DisplayName', 'Экспериментальные данные');

stairs(x0, f_nn, 'b', 'DisplayName', 'Метод ближайшего соседа');

plot(x0, f_linear, 'g', 'DisplayName', 'Линейная интерполяция');

plot(x0, f_lagrange, 'm', 'DisplayName', 'Полином Лагранжа');

title('Интерполяция: Ближайший сосед, Линейная, Лагранж');
legend;
hold off;

subplot(2, 1, 2);
hold on;
plot(x0, f0, 'k', 'LineWidth', 1.5, 'DisplayName', 'Аналитическая функция');

stem(x1, f1, 'r', 'filled', 'DisplayName', 'Экспериментальные данные');

plot(x0, f_newton_forward, 'b', 'DisplayName', 'Метод Ньютона вперед');

plot(x0, f_newton_backward, 'g', 'DisplayName', 'Метод Ньютона назад');

title('Интерполяция: Ньютон вперед, Ньютон назад');
legend;
hold off;


%1.3
function func = func(x)
    func = sin(2 * x);
end

%2.1
function f = sem_04_nn(x1, f1, x0)
    n = length(x0);
    f = zeros(1, n);
    for i = 1 : n
        [~, idx] = min(abs(x1 - x0(i)));
        f(i) = f1(idx);
    end
end

function f = sem_04_linear(x1, f1, x0)
    n = length(x0);
    f = zeros(1, n);
    for i = 1 : n
        if x0(i) <= x1(1)
            f(i) = f1(1);
        elseif x0(i) >= x1(end)
            f(i) = f1(end);
        else
            j = find(x1 <= x0(i), 1, 'last');
            f(i) = f1(j) + (f1(j + 1) - f1(j)) * (x0(i) - x1(j)) / (x1(j + 1) - x1(j));
        end
    end
end

function f = sem_04_lagrange(x1, f1, x0)
    n = length(x0);
    m = length(x1);
    f = zeros(1, n);
    for i = 1 : n
        L = ones(1, m);
        for j = 1 : m
            for k = 1 : m
                if k ~= j
                    L(j) = L(j) * (x0(i) - x1(k)) / (x1(j) - x1(k));
                end
            end
        end
        f(i) = sum(f1.* L);
    end
end

function f = sem_04_newton_forward(x1, f1, x0, n)
    m = length(x1);
    h = x1(2) - x1(1);

    diff_table = zeros(m, m);
    diff_table(:, 1) = f1';
    for j = 2 : m
        for i = 1 : (m - j + 1)
            diff_table(i, j) = diff_table(i + 1, j - 1) - diff_table(i, j - 1);
        end
    end
    
    f = zeros(1, length(x0));
    for i = 1 : length(x0)
        u = (x0(i) - x1(1)) / h;
        interp_val = diff_table(1, 1);
        u_term = 1;
        for j = 1 : n
            u_term = u_term * (u - j + 1) / j;
            interp_val = interp_val + u_term * diff_table(1, j + 1);
        end
        f(i) = interp_val;
    end
end

function f = sem_04_newton_backward(x1, f1, x0, n)
    m = length(x1);
    h = x1(2) - x1(1);
    
    diff_table = zeros(m, m);
    diff_table(:, 1) = f1';
    for j = 2 : m
        for i = 1 : (m - j + 1)
            diff_table(i, j) = diff_table(i + 1, j - 1) - diff_table(i, j - 1);
        end
    end
    
    f = zeros(1, length(x0));
    for i = 1 : length(x0)
        u = (x0(i) - x1(end)) / h;
        interp_val = diff_table(end, 1);
        u_term = 1;
        for j = 1 : n
            u_term = u_term * (u + j - 1) / j;
            interp_val = interp_val + u_term * diff_table(end - j, j + 1);
        end
        f(i) = interp_val;
    end
end
