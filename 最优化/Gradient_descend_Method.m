clear;
% close all;
clc

% 定义目标函数
f = @(x) x(1)^2 + x(2)^2 - 2*x(1)*x(2) + sin(x(1)) + cos(x(2));

% 定义目标函数的梯度
grad_f = @(x) [2*x(1) - 2*x(2) + cos(x(1)); 2*x(2) - 2*x(1) - sin(x(2))];




lr = 0.1;
max_iter = 100;

tol = 1e-1;
x0 = [20;-20];
history_x = zeros(2,max_iter);
history_f = zeros(1,max_iter);

x = x0;
for iteration = 1:max_iter
    gradient = grad_f(x);
    x_new = x - lr .* gradient;
    history_x(:,iteration) = x_new;
    history_f(iteration) = f(x_new);
    % 更新参数
    x = x_new;

    
    if norm(gradient) < tol
        break;
    end
end


% 可视化迭代过程
figure;
subplot(2, 1, 1);
plot(1:iteration, history_x(1, 1:iteration), '-o', 'LineWidth', 1.5);
hold on;
plot(1:iteration, history_x(2, 1:iteration), '-o', 'LineWidth', 1.5);
title('参数迭代过程');
legend('x(1)', 'x(2)');
xlabel('迭代次数');
ylabel('参数值');

subplot(2, 1, 2);
plot(1:iteration, history_f(1:iteration), '-o', 'LineWidth', 1.5);
title('目标函数值迭代过程');
xlabel('迭代次数');
ylabel('目标函数值');

% 显示最终结果
fprintf('最优解: x = [%f, %f]\n', x(1), x(2));
fprintf('f(x)的最优值: %f\n', f(x));
fprintf('迭代次数: %d\n', iteration);


% 目标函数可视化
axis_range = -20:1:20;
[X, Y] = meshgrid(axis_range, axis_range);
Z = arrayfun(@(x,y) f([x,y]), X, Y);
figure;
surf(X,Y,Z);
xlabel("X");
ylabel("Y");
zlabel("f");
grid on;
title("梯度下降法");

% 显示优化效果
hold on;
for j = 1:iteration
    
    plot3(history_x(1,j), history_x(2,j), f([history_x(1,j), history_x(2,j)]),'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
    % 连接到前一个点
    if j > 1
        line([history_x(1, j-1), history_x(1, j)], [history_x(2, j-1), history_x(2, j)], [f([history_x(1,j-1), history_x(2,j-1)]), f([history_x(1,j), history_x(2,j)])], 'Color', 'g');
    end
end


