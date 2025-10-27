clc; clear; close all;

% ============== 房间数据 ==============
% [emission, volume]
rooms = [
    1000, 100;
    5000, 300;
    15000, 300;
    90000, 60;
];

% 走廊 (emission=0, volume=500)
corridor = [0, 100];

% ============== 风量设定 ==============
q1 = 1; q2 = 2; q3 = 3;
q0 = q1 + q2 + q3;
assert(q3>q2 && q2>q1, 'q3 > q2 > q1 must hold');

% ============== 仿真设置 ==============
tspan = [0 1000];
initial_conditions = [0 0 0];
C0 = 0;
% ============== 阈值划分方案 ==============
threshold_schemes = [5000, 20000];

% ============== 走廊归属选项 ==============
corridor_options = {'0，Baseline','C2','C3'};

results = zeros(length(corridor_options), 10);
t_cell = cell(length(corridor_options),1);
C_cell = cell(length(corridor_options),1);

for s = 1:length(corridor_options)
    th1 = threshold_schemes(1);
    th2 = threshold_schemes(2);

    % 基于阈值的原始划分
    Q1 = sum(rooms(rooms(:,1) <= th1, 1));
    Q2 = sum(rooms(rooms(:,1) > th1 & rooms(:,1) <= th2, 1));
    Q3 = sum(rooms(rooms(:,1) > th2, 1));
    V1 = sum(rooms(rooms(:,1) <= th1, 2));
    V2 = sum(rooms(rooms(:,1) > th1 & rooms(:,1) <= th2, 2));
    V3 = sum(rooms(rooms(:,1) > th2, 2));

    % 加入走廊体积（除 baseline 外）
    if strcmp(corridor_options{s},'C2')
        V2 = V2 + corridor(2);
    elseif strcmp(corridor_options{s},'C3')
        V3 = V3 + corridor(2);
    end

    % 若某一区体积为0，跳过
    if any([V1,V2,V3]<=0)
        fprintf('Scheme %d (%s): infeasible (zero volume)\n', s, corridor_options{s});
        results(s,:) = [s, th1, th2, NaN, NaN, NaN, Inf, Inf, Inf, Inf];
        continue;
    end

    % ================== 集总稀释模型 ==================
    C1_eqn = @(t,C1) ( q0*(C0 - C1) + Q1 ) / V1;
    C2_eqn = @(t,C1,C2) ((q2+q3)*(C1 - C2) + Q2 ) / V2;
    C3_eqn = @(t,C2,C3) (q3 * (C2-C3)  + Q3)  / V3;

    % ODE 求解
    [t, C] = ode45(@(t,C)[C1_eqn(t,C(1)); C2_eqn(t,C(1),C(2)); C3_eqn(t,C(2),C(3))], ...
        tspan, initial_conditions);

    % 保存
    t_cell{s} = t;  C_cell{s} = C;

    % ================== 稳态及 90% 判据 ==================
    C1_inf = Q1 / q0;
    C2_inf = C1_inf + Q2 / (q2+q3);
    C3_inf = C2_inf + Q3 / q3;

    targets = 0.9 * [C1_inf, C2_inf, C3_inf];
    t90s = [local_time_to_reach(C(:,1),t,targets(1)), ...
            local_time_to_reach(C(:,2),t,targets(2)), ...
            local_time_to_reach(C(:,3),t,targets(3))];
    T_eq = max(t90s);

    results(s,:) = [s, th1, th2, C1_inf, C2_inf, C3_inf, t90s, T_eq];

    % ================== 绘图 ==================
    figure('Name', sprintf('Scheme %d: %s', s, corridor_options{s}));
    plot(t, C(:,1), 'r', 'LineWidth', 2); hold on;
    plot(t, C(:,2), 'g', 'LineWidth', 2);
    plot(t, C(:,3), 'b', 'LineWidth', 2);

    % 90%目标线
    yline(targets(1), 'r--', 'LineWidth', 2);
    yline(targets(2), 'g--', 'LineWidth', 2);
    yline(targets(3), 'b--', 'LineWidth', 2);

    % 标出最后达到90%的曲线点（黄色五角星 + 红边）
    if isfinite(T_eq)
        tol = 1e-6;
        idx_last = find(abs(t90s - T_eq) <= tol);
        for k = idx_last
            C_at_Teq = interp1(t, C(:,k), T_eq, 'linear', 'extrap');
            plot(T_eq, C_at_Teq, 'p', 'MarkerSize', 30, 'MarkerFaceColor', 'y', ...
                'MarkerEdgeColor', 'r', 'LineWidth', 1.8);
        end
        ymax = max(interp1(t, C, T_eq));
        text(T_eq, ymax*1.05, sprintf('T_{eq} = %.0f s', T_eq), ...
            'FontName','Arial','FontSize',30,'FontWeight','bold', ...
            'HorizontalAlignment','center');
    end

    % 格式统一
    xlabel('Time (s)', 'FontName','Arial', 'FontSize',25,'FontWeight','bold');
    ylabel('Concentration (Bq/m3)', 'FontName','Arial', 'FontSize',25,'FontWeight','bold');
    lg = legend({'C1','C2','C3','C1 90%','C2 90%','C3 90%'}, 'Location','best','FontSize',25);
    set(lg, 'FontName','Arial', 'FontSize',30,'FontWeight','bold');
    set(gca, 'FontName','Arial', 'FontSize',30,'FontWeight','bold','LineWidth', 1.5);
    grid on; box on;
end

% ============== 输出结果表格 ==============
T = array2table(results, ...
    'VariableNames', {'Scheme','Th1','Th2','C1_inf','C2_inf','C3_inf','t90_C1','t90_C2','t90_C3','T_eq'});
T.SchemeName = corridor_options';
disp(T);

% 找出最佳方案
T_eqs = results(:,10);
finite_idx = isfinite(T_eqs);
if any(finite_idx)
    [best_Teq, relidx] = min(T_eqs(finite_idx));
    all_idx = find(finite_idx);
    best_row_idx = all_idx(relidx);
    fprintf('\n>>>> Best scheme: %s, T_eq = %.2f s\n', corridor_options{best_row_idx}, best_Teq);
end

% ===================== 辅助函数 =====================
function tt = local_time_to_reach(Cvec, tvec, target)
    if isnan(target) || isinf(target), tt = Inf; return; end
    if target <= 0, tt = tvec(1); return; end
    idx = find(Cvec >= target, 1, 'first');
    if isempty(idx)
        tt = Inf;
    elseif idx == 1
        tt = tvec(1);
    else
        Clo = Cvec(idx-1); Chi = Cvec(idx);
        tlo = tvec(idx-1); thi = tvec(idx);
        tt = tlo + (target - Clo)*(thi-tlo)/(Chi-Clo);
    end
end
QQ= V1*C1_inf+C2_inf*V2+C3_inf*V3;
disp(QQ);