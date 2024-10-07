
%3.1 3.2
A = randn(5);
X1 = expm(A);
X2 = sem_01_expm(A, 5);
max_val = max(max(abs(X1 - X2)));
% fprintf('Масимальное расхождение: %d\n', max_val);


% 3.3
analyse_func_01_expm_results()


% 4.1
med = analyse_func_01_expm_calc_time();
fprintf("mid_time: %d\n", med);


function result = analyse_func_01_expm_calc_time()
    time_result = zeros(5000, 1);
    for i=1:5000
        A = randn(5);
        X1 = expm(A);
        tic;
        X2 = sem_01_expm(A, 8);
        time_result(i) = toc;
    end
    figure Name 'time_data' ;   
    plot(1:4000, time_result(501:1:4500));
    result = median(time_result(501:1:4500));
end


function analyse_func_01_expm_results()
    y = zeros(15, 1);
    A = randn(5);
    for k=1:15
        X1 = expm(A);
        X2 = sem_01_expm(A, k);
        
        max_val = max(max(abs(X1 - X2)));
        y(k) = max_val;
        % fprintf('Масимальное расхождение: %d\n', max_val);
    end
    figure Name 'accuracy_graph';
    bar(1:15, y);
end



function result = sem_01_expm(A, k)
    assert (k >= 1);
    result = zeros(size(A));
    for i=0:k
        result = result + 1/factorial(i) * A^i;
    end
end