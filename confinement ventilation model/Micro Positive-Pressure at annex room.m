clc; clear; close all;

% ============== 房间数据 ==============
% [emission (m_dot), volume (m³)]
rooms = [
   1000, 100;
    5000, 300;
    15000, 300;
    90000, 60;
];
corridor = [0]; % 走廊 (emission=0, volume=500)

% ============== 风量设定 ==============
q1 = 1; q2 = 2; q3 = 3;
q0 = 1.5; % 环境渗透量
qx1 = 0.5; qx2 = 1.5; % 渗透风量
assert(q3>q2 && q2>q1, '风量关系错误: q3>q2>q1');

% ============== 仿真设置 ==============
tspan = [0 1000];
initial_conditions = [0 0 0];
C0 = 0;

% ============== 阈值划分方案 ==============
th1 = 5000; th2 = 20000;

% ============== 两个方案：不含走廊 vs 含走廊 ==============
schemes = {'Baseline (no corridor)','With corridor'};
results = zeros(length(schemes), 10);

for s = 1:length(schemes)
    % 划分区域
    Q1 = sum(rooms(rooms(:,1) <= th1, 1));
    Q2 = sum(rooms(rooms(:,1) > th1 & rooms(:,1) <= th2, 1));
    Q3 = sum(rooms(rooms(:,1) > th2, 1));
    V1 = sum(rooms(rooms(:,1) <= th1, 2));
    V2 = sum(rooms(rooms(:,1) > th1 & rooms(:,1) <= th2, 2));
    V3 = sum(rooms(rooms(:,1) > th2, 2));
 

    % ============= 新 ODE 模型 =============
    dC1 = @(t,C1)( q0*(C0 - C1) + Q1 ) / V1;
    dC2 = @(t,C1,C2) ( (q2 - qx2)*(C1 - C2) + qx2*(C0 - C2) + Q2 ) / V2;
    dC3 = @(t, C3) ( q3*(C0 - C3) + Q3 ) / V3;

    [t, C] = ode45(@(t,C)[dC1(t,C(1));
                           dC2(t,C(1),C(2));
                           dC3(t,C(3))], tspan, initial_conditions);

    % 稳态浓度（t→∞）
    C1_inf = C0 + Q1/q0;
    C2_inf = ((q2 - qx2)*C1_inf + qx2*C0 + Q2)/q2;
    C3_inf = C0 + Q3/q3;
    targets = 0.9*[C1_inf, C2_inf, C3_inf];

    t90s = [local_time_to_reach(C(:,1),t,targets(1)), ...
            local_time_to_reach(C(:,2),t,targets(2)), ...
            local_time_to_reach(C(:,3),t,targets(3))];
    T_eq = max(t90s);

    results(s,:) = [s, th1, th2, C1_inf, C2_inf, C3_inf, t90s, T_eq];
   
    % ============= 绘图 =============
    figure('Name', schemes{s});
    plot(t,C(:,1),'r','LineWidth',2); hold on;
    plot(t,C(:,2),'g','LineWidth',2);
    plot(t,C(:,3),'b','LineWidth',2);
  yline(targets(1), 'r--', 'LineWidth', 2);
    yline(targets(2), 'g--', 'LineWidth', 2);
    yline(targets(3), 'b--', 'LineWidth', 2);

    % 黄色五角星标注最后平衡区
    if isfinite(T_eq)
        idx_last = find(abs(t90s - T_eq) < 1e-6);
        for k = idx_last
            C_at = interp1(t,C(:,k),T_eq);
            plot(T_eq, C_at, 'p', 'MarkerSize', 30, 'MarkerFaceColor', 'y', ...
                 'MarkerEdgeColor', 'r', 'LineWidth', 1.8);
        end
        ymax = max(interp1(t,C,T_eq));
        text(T_eq, ymax*1.05, sprintf('T_{eq}=%.0f s',T_eq), ...
            'FontName','Arial','FontSize',30,'FontWeight','bold', ...
            'HorizontalAlignment','center');
    end

    xlabel('Time (s)', 'FontName','Arial', 'FontSize',25,'FontWeight','bold');
    ylabel('Concentration (Bq/m3)', 'FontName','Arial', 'FontSize',25,'FontWeight','bold');
    lg = legend({'C1','C2','C3','C1 90%','C2 90%','C3 90%'}, 'Location','best','FontSize',25);
    set(lg, 'FontName','Arial', 'FontSize',30,'FontWeight','bold');
    set(gca, 'FontName','Arial', 'FontSize',30,'FontWeight','bold','LineWidth', 1.5);
    grid on; box on;
end

% ============= 输出结果表格 ==============
T = array2table(results, ...
    'VariableNames',{'Scheme','Th1','Th2','C1_inf','C2_inf','C3_inf','t90_C1','t90_C2','t90_C3','T_eq'});
T.SchemeName = schemes';
disp(T);

% ===================== 辅助函数 =====================
function tt = local_time_to_reach(Cvec, tvec, target)
    if isnan(target)||isinf(target), tt=Inf; return; end
    idx = find(Cvec>=target,1,'first');
    if isempty(idx)
        tt = Inf;
    elseif idx==1
        tt = tvec(1);
    else
        Clo=Cvec(idx-1); Chi=Cvec(idx);
        tlo=tvec(idx-1); thi=tvec(idx);
        tt=tlo+(target-Clo)*(thi-tlo)/(Chi-Clo);
    end
end
QQ= V1*C1_inf+C2_inf*V2+C3_inf*V3;
disp(QQ);