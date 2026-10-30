%%
% 作者：雷宇杰 | 班级：储能231 | 学号：206231313 | 本人任务：模型A，节点4/9/14
%% 第一周作业：计及配电网二阶锥潮流与充放电互斥约束 (MISOCP) 的储能日前优化调度
% 课程名称：《储能优化配置与调度技术》课程设计
%% 0. 模板使用说明与核心数学模型
% *使用说明*：
%%
% * 请将本模板复制到 |ESS-Course-Submissions| 仓库下你的 【学号_姓名拼音】 文件夹中；
% * 重命名为 【姓名拼音_学号_Week1_Task3.mlx】，并修改第 2 节中的 |Stu_ID| 为本人真实 9 位学号；
% * 仔细阅读各节的 LaTeX 数学模型说明，在标有 |% 【已完成】| 的代码块处补全核心建模语句，并导出对应的 |.m| 纯文本脚本一并提交。
%%
% *核心数学模型*：
%
% *1) 储能系统荷电状态（SoC）动态递推方程与充放电互斥约束*（已知单站额定功率 $P_i^{\mathrm{rated}}$、额定容量 $E_i^{\mathrm{rated}}$
% 为常数，引入二进制互斥变量 $u_{i,t}^{\mathrm{ch}} \in \{0,1\}$）：
%
% $$\mathrm{SoC}_{i,t} = \mathrm{SoC}_{i,t-1} + \frac{\eta_{\mathrm{ch}} P_{i,t}^{\mathrm{ch}}
% \Delta t}{E_i^{\mathrm{rated}}} - \frac{P_{i,t}^{\mathrm{dis}} \Delta t}{\eta_{\mathrm{dis}}
% E_i^{\mathrm{rated}}}, \quad \mathrm{SoC}_{\min} \le \mathrm{SoC}_{i,t} \le
% \mathrm{SoC}_{\max}$$
%
% $$0 \le P_{i,t}^{\mathrm{ch}} \le u_{i,t}^{\mathrm{ch}} P_i^{\mathrm{rated}},
% \quad 0 \le P_{i,t}^{\mathrm{dis}} \le \left(1 - u_{i,t}^{\mathrm{ch}}\right)
% P_i^{\mathrm{rated}}, \quad u_{i,t}^{\mathrm{ch}} \in \{0,1\}$$
%
% *2) 配电网 DistFlow 支路潮流与旋转二阶锥松弛方程*（令 $U_{i,t} = V_{i,t}^2, L_{ij,t} = I_{ij,t}^2$）：
%
% $$P_{ij,t} - r_{ij} L_{ij,t} - \sum_{k\in\delta(j)} P_{jk,t} = P_{j,t}^{\mathrm{load}}
% - P_{j,t}^{\mathrm{pv}} + P_{j,t}^{\mathrm{ch}} - P_{j,t}^{\mathrm{dis}}$$
%
% $$Q_{ij,t} - x_{ij} L_{ij,t} - \sum_{k\in\delta(j)} Q_{jk,t} = Q_{j,t}^{\mathrm{load}}$$
%
% $$U_{j,t} = U_{i,t} - 2\left(r_{ij} P_{ij,t} + x_{ij} Q_{ij,t}\right) + \left(r_{ij}^2
% + x_{ij}^2\right) L_{ij,t}$$
%
% $$\sqrt{\left(2P_{ij,t}\right)^2 + \left(2Q_{ij,t}\right)^2 + \left(L_{ij,t}
% - U_{i,t}\right)^2} \le L_{ij,t} + U_{i,t}, \quad \left(V_{\min}\right)^2 \le
% U_{i,t} \le \left(V_{\max}\right)^2$$

clear; clc; close all;
objective_number=0; % 0按学号分配；1/2/3分别选择A/B/C，保持算例与节点不变。
% 雷宇杰，储能231，206231313。Live Editor可能在临时目录求值，须核对真实源文件。
source_stem='Leiyujie_206231313_Week1_Task3';
script_dir=fileparts(mfilename('fullpath'));
if isempty(script_dir),script_dir=pwd;end
if ~isfile(fullfile(script_dir,[source_stem '.mlx']))
    candidates={getenv('ESS_WEEK1_DIR'),pwd,fileparts(which([source_stem '.m'])),...
        fileparts(which([source_stem '.mlx']))};
    active_file=matlab.desktop.editor.getActiveFilename;
    [active_dir,active_stem]=fileparts(active_file);
    if strcmp(active_stem,source_stem),candidates{end+1}=active_dir;end
    ancestor=pwd;
    while true
        candidates{end+1}=fullfile(ancestor,'ESS-Course-Submissions','206231313_Leiyujie');
        parent=fileparts(ancestor);if strcmp(parent,ancestor),break;end
        ancestor=parent;
    end
    script_dir='';
    for candidate_index=1:numel(candidates)
        candidate=candidates{candidate_index};
        if ~isempty(candidate) && isfile(fullfile(candidate,[source_stem '.mlx']))
            script_dir=candidate;break;
        end
    end
end
assert(~isempty(script_dir) && isfile(fullfile(script_dir,[source_stem '.mlx'])),...
    '未定位真实作业目录。请把MATLAB当前文件夹切换到个人作业目录再运行。');
project_root=fileparts(fileparts(script_dir));
if isfile(fullfile(project_root,'setup_env.m'))
    addpath(project_root);setup_env;
else
    course_dir=fullfile(project_root,'ESS-Course-Materials');
    if ~isfolder(course_dir) && exist('case_nongcun','file')~=0
        course_dir=fileparts(fileparts(which('case_nongcun')));
    end
    if isfolder(fullfile(course_dir,'cases')),addpath(fullfile(course_dir,'cases'));end
    assert(exist('runpf','file')~=0 && exist('sdpvar','file')~=0,...
        '请先配置MATPOWER、YALMIP及MOSEK/Gurobi。');
end
addpath(script_dir);
cd(script_dir);
yalmip('clear');
fprintf('OUTPUT_DIRECTORY=%s\n',script_dir);

%% 1. 自动检测算例库是否已在 MATLAB 搜索路径中

if ~exist('case_nongcun.m', 'file')
    % 若不在路径中，弹窗提示学生选择一次本地的 ESS-Course-Materials/cases 文件夹
    case_dir = uigetdir(pwd, '请选择教学资源库中的 cases 文件夹路径');
    addpath(case_dir);
end
%% 2. 学号映射与网络参数读取
% *【学生修改区】*：请将下方 |Stu_ID| 修改为你本人的 9 位学号！程序将根据末两位 $\mathrm{XY}$ 自动计算负荷缩放乘子 $\alpha
% = 1 + 0.002 \times \mathrm{XY}$ 并分配专属算例。

Stu_ID  = 206231313;
XY      = mod(Stu_ID, 100);
X_digit = floor(XY / 10);
Y_digit = mod(XY, 10);
alpha   = 1 + 0.002 * XY;

if ismember(Y_digit, [0, 1, 2, 3])
    case_name = 'case_nongcun';
elseif ismember(Y_digit, [4, 5, 6])
    case_name = 'case_chengshi_A';
else
    case_name = 'case_kaifaqu_A';
end
mpc = loadcase(case_name);
load('profiles_24h.mat');

if ismember(X_digit, [0, 1, 2])
    ess_buses  = [4, 9, 14];
    model_type = 'A'; % 模型 A：最大化新能源消纳（含微网损正则项）
elseif ismember(X_digit, [3, 4, 5])
    ess_buses  = [5, 11, 16];
    model_type = 'B'; % 模型 B：最小化全网有功网损
else
    ess_buses  = [6, 13, 17];
    model_type = 'C'; % 模型 C：最小化主网购电成本
end
pv_buses = ess_buses;
if objective_number~=0
 assert(ismember(objective_number,1:3));
 objective_names='ABC';model_type=objective_names(objective_number);
end
%% 3. 标幺值基准转换与设备参数设置

Sbase = mpc.baseMVA;
N_bus = size(mpc.bus, 1);
N_br  = size(mpc.branch, 1);
N_ess = length(ess_buses);
N_pv  = length(pv_buses);

fbus = mpc.branch(:, 1);
tbus = mpc.branch(:, 2);
r_br = mpc.branch(:, 3);
x_br = mpc.branch(:, 4);

Pd_pu       = (mpc.bus(:, 3) * alpha / Sbase) * load_curve_pu';
Qd_pu       = (mpc.bus(:, 4) * alpha / Sbase) * load_curve_pu';
P_pv_max_pu = (1.2 / Sbase) * repmat(pv_curve_pu', N_pv, 1);

P_ess_rated_pu = 1.5 / Sbase; % 单站固定额定功率 1.5 MW (标幺值)
E_ess_rated_pu = 3.0 / Sbase; % 单站固定额定容量 3.0 MWh (标幺值)
eta_ch  = 0.90;
eta_dis = 0.95;
SoC_min = 0.20;
SoC_max = 0.80;
SoC_0   = 0.50;
%% 4. 声明 YALMIP 决策变量

U_bus  = sdpvar(N_bus, T, 'full'); % 节点电压平方 U_{i,t} = V_{i,t}^2
L_br   = sdpvar(N_br,  T, 'full'); % 支路电流平方 L_{ij,t} = I_{ij,t}^2
P_br   = sdpvar(N_br,  T, 'full'); % 支路首端有功功率 (p.u.)
Q_br   = sdpvar(N_br,  T, 'full'); % 支路首端无功功率 (p.u.)
P_grid = sdpvar(1,     T, 'full'); % 关口有功功率 (p.u.)
Q_grid = sdpvar(1,     T, 'full'); % 关口无功功率 (p.u.)
P_ch   = sdpvar(N_ess, T, 'full'); % 储能充电有功功率 (p.u.)
P_dis  = sdpvar(N_ess, T, 'full'); % 储能放电有功功率 (p.u.)
u_ch   = binvar(N_ess, T, 'full'); % 储能充放电互斥二进制变量 u_{i,t}^{ch} \in {0,1}
SoC    = sdpvar(N_ess, T, 'full'); % 储能荷电状态 SoC_{i,t}
P_pv   = sdpvar(N_pv,  T, 'full'); % 分布式光伏实际出力 (p.u.)
%% 5. 构建约束条件
%% 5.1 储能充放电互斥与分布式光伏运行约束

Constraints = [];
Constraints = [Constraints, 0 <= P_ch  <= P_ess_rated_pu * u_ch];
Constraints = [Constraints, 0 <= P_dis <= P_ess_rated_pu * (1 - u_ch)];
Constraints = [Constraints, SoC_min <= SoC <= SoC_max];
Constraints = [Constraints, SoC(:, T) == SoC_0];
Constraints = [Constraints, 0 <= P_pv <= P_pv_max_pu];

for t = 1:T
    if t == 1
        soc_prev = repmat(SoC_0, N_ess, 1);
    else
        soc_prev = SoC(:, t - 1);
    end
    % 【完成 1】储能递推：储能荷电状态 SoC(:, t) 的跨时段耦合递推方程：
    Constraints = [Constraints, SoC(:, t) == soc_prev + ...
        (eta_ch * P_ch(:, t) - P_dis(:, t) / eta_dis) / E_ess_rated_pu];
end
%% 5.2 节点功率平衡与 DistFlow 二阶锥潮流约束

A_from  = sparse(fbus, 1:N_br, 1, N_bus, N_br);
A_to    = sparse(tbus, 1:N_br, 1, N_bus, N_br);
A_ess   = sparse(ess_buses, 1:N_ess, 1, N_bus, N_ess);
A_pv    = sparse(pv_buses,  1:N_pv,  1, N_bus, N_pv);
A_slack = sparse(1, 1, 1, N_bus, 1);

P_loss_br = diag(r_br) * L_br;
Q_loss_br = diag(x_br) * L_br;

Constraints = [Constraints, A_to * (P_br - P_loss_br) - A_from * P_br + ...
    A_slack * P_grid + A_pv * P_pv + A_ess * (P_dis - P_ch) == Pd_pu];
Constraints = [Constraints, A_to * (Q_br - Q_loss_br) - A_from * Q_br + ...
    A_slack * Q_grid == Qd_pu];

Constraints = [Constraints, U_bus(1, :) == 1.0];
Constraints = [Constraints, (0.95^2) <= U_bus(2:end, :) <= (1.05^2)];
Constraints = [Constraints, L_br >= 0];

for b = 1:N_br
    i = fbus(b);
    j = tbus(b);
    % 【完成 2】DistFlow： DistFlow 支路电压平方降落方程与 cone 二阶锥松弛约束：
    Constraints = [Constraints, U_bus(j, :) == U_bus(i, :) - ...
        2 * (r_br(b) * P_br(b, :) + x_br(b) * Q_br(b, :)) + ...
        (r_br(b)^2 + x_br(b)^2) * L_br(b, :)];
    for t = 1:T
        Constraints = [Constraints, cone([2 * P_br(b, t); 2 * Q_br(b, t); ...
            L_br(b, t) - U_bus(i, t)], L_br(b, t) + U_bus(i, t))];
    end
end
%% 6. 构建目标函数并调用 MISOCP 求解器 (MOSEK / Gurobi)
% 本模型同时包含充放电互斥二进制变量 $u_{i,t}^{\mathrm{ch}}$ 与支路二阶锥约束，属于混合整数二阶锥规划（MISOCP），需调用
% MOSEK 或 Gurobi 求解。若使用 CVX 随附的 MOSEK，请仅将 |cvx/mosek| 子文件夹加入 MATLAB 路径（|addpath(genpath('.../cvx/mosek'))|），切勿运行
% |cvx_setup| 或添加 |cvx| 根目录。

Total_Loss_MWh    = Sbase * sum(sum(P_loss_br));
Total_Curtail_MWh = Sbase * sum(sum(P_pv_max_pu - P_pv));
Total_Cost_Yuan   = sum((Sbase * P_grid') .* price_MWh);

% 【完成 3】按本人学号映射的 model_type ('A' / 'B' / 'C') 选择目标函数：
if strcmp(model_type, 'A')
    Objective = Total_Curtail_MWh + 1e-2 * Total_Loss_MWh;
elseif strcmp(model_type, 'B')
    Objective = Total_Loss_MWh + 1e-4 * sum(sum(P_ch + P_dis));
else
    Objective = Total_Cost_Yuan;
end

solver_candidates={};
if exist('mosekopt','file')~=0,solver_candidates{end+1}='mosek';end
if exist('gurobi','file')~=0,solver_candidates{end+1}='gurobi';end
assert(~isempty(solver_candidates),...
    '本任务含充放电二进制变量和二阶锥约束，需配置 MOSEK 或 Gurobi。');
messages=cell(0,1);solver_succeeded=false;tic;
for k=1:numel(solver_candidates)
    solver_name=solver_candidates{k};
    ops=sdpsettings('solver',solver_name,'verbose',1,'savesolveroutput',1);
    % 由本函数捕获原始异常，既保留错误原因，也继续尝试下一个求解器。
    ops.debug=1;
    if strcmp(solver_name,'mosek')
        ops.mosek.MSK_DPAR_MIO_MAX_TIME=600;
        ops.mosek.MSK_DPAR_MIO_TOL_REL_GAP=1e-6;
        ops.mosek.MSK_DPAR_MIO_TOL_ABS_GAP=1e-8;
        ops.mosek.MSK_DPAR_MIO_TOL_FEAS=1e-9;
        ops.mosek.MSK_DPAR_MIO_TOL_ABS_RELAX_INT=1e-9;
        ops.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-10;
        ops.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-10;
        ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-10;
    else
        % 旧版 YALMIP 默认值为 -1，新版 Gurobi 要求非负。
        ops.gurobi.TuneTimeLimit=0;
        ops.gurobi.TimeLimit=120;
        ops.gurobi.MIPGap=1e-6;
        ops.gurobi.FeasibilityTol=1e-9;
        ops.gurobi.IntFeasTol=1e-9;
        ops.gurobi.OptimalityTol=1e-8;
        ops.gurobi.BarQCPConvTol=1e-9;
    end
    try
        sol=optimize(Constraints,Objective,ops);
    catch solver_error
        messages{end+1}=sprintf('%s [%s]: %s',solver_name,...
            solver_error.identifier,solver_error.message); %#ok<AGROW>
        continue
    end
    if sol.problem==0,solve_seconds=toc;solver_succeeded=true;break;end
    detail=sol.info;
    if isfield(sol,'solveroutput') && ~isempty(sol.solveroutput)
        % MOSEK 的 rcode/rmsg/rcodestr 及终止状态保留在原始输出中。
        detail=sprintf('%s\n原始求解器输出：\n%s',detail,evalc('disp(sol.solveroutput)'));
        if strcmp(solver_name,'mosek') && isfield(sol.solveroutput,'res')
            raw=sol.solveroutput.res;
            detail=sprintf('%s\nMOSEK 返回信息：\n%s',detail,evalc('disp(raw)'));
            if isfield(raw,'sol') && isfield(raw.sol,'int')
                integer_solution=raw.sol.int;
                if isfield(integer_solution,'prosta')
                    detail=sprintf('%s\n整数问题状态：%s',detail,integer_solution.prosta);
                end
                if isfield(integer_solution,'solsta')
                    detail=sprintf('%s\n整数解状态：%s',detail,integer_solution.solsta);
                end
            end
            if isfield(raw,'rcode') && raw.rcode==100001
                detail=sprintf('%s\nMOSEK 已达到本次设置的 %g 秒时间上限。',...
                    detail,ops.mosek.MSK_DPAR_MIO_MAX_TIME);
            end
        end
    end
    messages{end+1}=[solver_name ': ' detail]; %#ok<AGROW>
end
if ~solver_succeeded
    error('MISOCP 求解失败：%s',strjoin(messages,'; '));
end
%% 7. 二阶锥松弛精确性分析、交流潮流回代与《中国电机工程学报》规范绘图
% 请注意图片字体与符号规范：中文采用宋体，西文与数字采用 Times New Roman，变量符号采用斜体，单位与下标说明采用正体。

socp_gap     = max(max(abs(value(L_br) .* value(U_bus(fbus, :)) - (value(P_br).^2 + value(Q_br).^2))));
simul_err_MW = max(max(min(value(P_ch), value(P_dis)))) * Sbase;

V_socp_pu     = sqrt(value(U_bus));
V_ac_ess      = zeros(N_bus, T);
Ploss_ac_MWh  = 0;
P_net_ess_MW  = Sbase * value(P_dis - P_ch);
P_pv_MW_val   = Sbase * value(P_pv);
mpopt = mpoption('verbose', 0, 'out.all', 0);
for t = 1:T
    mpc_t = mpc;
    mpc_t.bus(:, 3) = Pd_pu(:, t) * Sbase;
    mpc_t.bus(:, 4) = Qd_pu(:, t) * Sbase;
    mpc_t.bus(pv_buses,  3) = mpc_t.bus(pv_buses,  3) - P_pv_MW_val(:, t);
    mpc_t.bus(ess_buses, 3) = mpc_t.bus(ess_buses, 3) - P_net_ess_MW(:, t);
    res_t = runpf(mpc_t, mpopt);
    assert(res_t.success==1,'第%d时段交流潮流未收敛。',t);
    V_ac_ess(:, t) = res_t.bus(:, 8);
    Ploss_ac_MWh   = Ploss_ac_MWh + sum(real(get_losses(res_t)));
end
max_ac_v_err = max(max(abs(V_socp_pu - V_ac_ess)));

fprintf('\n================ 第一周：配电网储能 MISOCP 优化调度结果 ================\n');
fprintf('学号: %d | 算例: %s | 模型: %s | 求解器: %s\n', Stu_ID, case_name, model_type, solver_name);
fprintf('模型全天电能损耗: %.4f MWh | 交流潮流电能损耗: %.4f MWh | 弃光量: %.4f MWh | 关口购电成本: %.2f 元\n', ...
    value(Total_Loss_MWh), Ploss_ac_MWh, value(Total_Curtail_MWh), value(Total_Cost_Yuan));
fprintf('同时充放电最大重叠功率: %.2e MW | 最大二阶锥松弛间隙: %.2e p.u.\n', simul_err_MW, socp_gap);
fprintf('交流潮流节点电压范围: [%.4f, %.4f] p.u. | 模型与交流电压最大偏差: %.2e p.u.\n', ...
    min(V_ac_ess(:)), max(V_ac_ess(:)), max_ac_v_err);
if socp_gap > 1e-4
    fprintf('[教学提示] 当前模型(%s)目标函数不含网损正惩罚，高光伏时段节点电压触及上限，出现二阶锥松弛不精确与交流电压越限特例。\n', model_type);
end
fprintf('========================================================================\n');

P_net_MW_val = Sbase * value(P_dis - P_ch)';
SoC_pct_val  = value(SoC)' * 100;
line_styles  = {'-o', '-s', '-^'};
colors       = [0.00, 0.45, 0.74; 0.85, 0.33, 0.10; 0.47, 0.67, 0.19];

fig1 = figure('Name', '配电网储能日前SOCP调度结果', 'Color', 'w', 'Position', [120, 100, 720, 500]);
if isprop(fig1, 'Theme'), theme(fig1, 'light'); end

subplot(2, 1, 1);
hold on;
for k = 1:N_ess
    plot(1:T, P_net_MW_val(:, k), line_styles{k}, 'Color', colors(k, :), ...
        'LineWidth', 1.6, 'MarkerSize', 5, 'MarkerFaceColor', colors(k, :));
end
yline(0, '--k', 'LineWidth', 0.9);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontName', 'Times New Roman', 'FontSize', 10.5, 'LineWidth', 1.0, ...
    'Box', 'on', 'TickDir', 'in', 'XTick', 0:4:24, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'GridColor', [0.4, 0.4, 0.4]);
xlim([1, T]);
ylim([min([P_net_MW_val(:); -0.5]) - 0.25, max([P_net_MW_val(:); 0.5]) + 0.55]);
xlabel('\fontname{宋体}时间 \fontname{Times New Roman}{\itt}\rm{/h}', 'FontSize', 11.5, 'Color', 'k', 'Interpreter', 'tex');
ylabel('\fontname{宋体}储能净放电功率 \fontname{Times New Roman}{\itP}_{\rm{ess}}\rm{/MW}', 'FontSize', 11.5, 'Color', 'k', 'Interpreter', 'tex');
title(sprintf('\\fontname{宋体}(a) 各候选节点储能日前调度有功出力 (模型 \\fontname{Times New Roman}%s)', model_type), ...
    'FontSize', 11.5, 'FontWeight', 'normal', 'Color', 'k', 'Interpreter', 'tex');
leg_str_a = arrayfun(@(b) sprintf('\\fontname{宋体}节点 \\fontname{Times New Roman}%d \\fontname{宋体}储能', b), ...
    ess_buses, 'UniformOutput', false);
leg1 = legend(leg_str_a, 'Location', 'north', 'Orientation', 'horizontal', 'FontSize', 10, 'Interpreter', 'tex', 'Box', 'on');
set(leg1, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
grid on;

subplot(2, 1, 2);
hold on;
for k = 1:N_ess
    plot(1:T, SoC_pct_val(:, k), line_styles{k}, 'Color', colors(k, :), ...
        'LineWidth', 1.6, 'MarkerSize', 5, 'MarkerFaceColor', colors(k, :));
end
yline(SoC_min * 100, ':k', 'LineWidth', 1.1);
yline(SoC_max * 100, ':r', 'LineWidth', 1.1);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontName', 'Times New Roman', 'FontSize', 10.5, 'LineWidth', 1.0, ...
    'Box', 'on', 'TickDir', 'in', 'XTick', 0:4:24, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'GridColor', [0.4, 0.4, 0.4]);
xlim([1, T]);
ylim([10, 98]);
xlabel('\fontname{宋体}时间 \fontname{Times New Roman}{\itt}\rm{/h}', 'FontSize', 11.5, 'Color', 'k', 'Interpreter', 'tex');
ylabel('\fontname{宋体}荷电状态 \fontname{Times New Roman}\rm{SoC}/%', 'FontSize', 11.5, 'Color', 'k', 'Interpreter', 'tex');
title('\fontname{宋体}(b) 各候选节点储能荷电状态 \fontname{Times New Roman}\rm{SoC} \fontname{宋体}时序变化曲线', ...
    'FontSize', 11.5, 'FontWeight', 'normal', 'Color', 'k', 'Interpreter', 'tex');
leg2 = legend(leg_str_a, 'Location', 'north', 'Orientation', 'horizontal', 'FontSize', 10, 'Interpreter', 'tex', 'Box', 'on');
set(leg2, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
grid on;

% 默认自动导出高质量白底 PNG (600 dpi)、矢量 PDF 与 MATLAB 源文件 FIG
%% 8. 数值审计、整数最优性证据与第一周结果归档
% 对同一组出力重新检查能量守恒、互斥、支路平衡、锥间隙与24h交流电压。
Dispatch_Data=struct('mpc',mpc,'Stu_ID',Stu_ID,'case_name',case_name,...
    'model_type',model_type,'solver',solver_name,'solve_seconds',solve_seconds,...
    'Sbase',Sbase,'dt',1,'ess_buses',ess_buses,'pv_buses',pv_buses,...
    'PratedMW',P_ess_rated_pu*Sbase,'EratedMWh',E_ess_rated_pu*Sbase,...
    'SoC0',SoC_0,'SoCmin',SoC_min,'SoCmax',SoC_max,...
    'eta_ch',eta_ch,'eta_dis',eta_dis,'price_MWh',price_MWh,...
    'PloadMW',Pd_pu*Sbase,'QloadMvar',Qd_pu*Sbase,...
    'PpvMW',Sbase*value(P_pv),'PpvMaxMW',P_pv_max_pu*Sbase,...
    'PchMW',Sbase*value(P_ch),'PdisMW',Sbase*value(P_dis),...
    'SoC',value(SoC),'u',value(u_ch),'U',value(U_bus),...
    'L',value(L_br),'P',value(P_br),'Q',value(Q_br),...
    'Pgrid',value(P_grid),'Qgrid',value(Q_grid));
[Dispatch_Check,V_dispatch_ac]=week1_check_dispatch(Dispatch_Data);
assert(Dispatch_Check.StorageConstraintsPassed && Dispatch_Check.LinearConstraintsPassed,...
    '储能或配电网约束检查失败。');
assert(Dispatch_Check.TeachingRequirementMet,'锥间隙或交流电压检查失败。');
if strcmp(model_type,'A'),assert(Dispatch_Check.Passed);end

IntegerStatus='UNAVAILABLE'; AbsoluteGap=NaN; RelativeGap=NaN; LowerBound=NaN;
ProvenOptimal=false; SolverRaw=sol.solveroutput;
if strcmp(solver_name,'mosek') && isfield(SolverRaw,'res')
    raw=SolverRaw.res;
    if isfield(raw,'sol') && isfield(raw.sol,'int')
        IntegerStatus=raw.sol.int.solsta;
    end
    if isfield(raw,'info') && isfield(raw.info,'MSK_DINF_MIO_OBJ_BOUND')
        % 求解器目标可能移除了常数；目标上下界之差不受常数影响。
        AbsoluteGap=abs(raw.info.MSK_DINF_MIO_OBJ_INT-raw.info.MSK_DINF_MIO_OBJ_BOUND);
        RelativeGap=AbsoluteGap/max(abs(value(Objective)),1e-12);
        LowerBound=value(Objective)-AbsoluteGap;
    end
    ProvenOptimal=strcmp(IntegerStatus,'INTEGER_OPTIMAL') && ...
        AbsoluteGap<=max(1e-8,1e-6*abs(value(Objective)));
end
Result_Week1=struct('StuID',Stu_ID,'Name','雷宇杰','NamePinyin','Leiyujie',...
    'Class','储能231','Stage','Week1','Case',case_name,'Alpha',alpha,...
    'Model',model_type,'ESS_Buses',ess_buses,'Solver',solver_name,...
    'MATLABVersion',version,'MaterialCommit','1d34d745ee234be0d27d8790dc95d687538fab2b',...
    'ObjectiveValue',value(Objective),'LineLoss_MWh',Dispatch_Check.LineLoss_MWh,...
    'ACLineLoss_MWh',Dispatch_Check.ACLineLoss_MWh,...
    'Curtailment_MWh',Dispatch_Check.Curtailment_MWh,'GridCost_Yuan',Dispatch_Check.GridCost_Yuan,...
    'MaxConeGap',Dispatch_Check.MaxConeGap,'Max_V_Err',Dispatch_Check.MaxACVoltageError,...
    'MinACVoltage',Dispatch_Check.MinACVoltage,'MaxACVoltage',Dispatch_Check.MaxACVoltage,...
    'SolveSeconds',solve_seconds,'IntegerStatus',IntegerStatus,...
    'AbsoluteGap',AbsoluteGap,'RelativeGap',RelativeGap,'LowerBound',LowerBound,...
    'ProvenOptimal',ProvenOptimal,'FeasibilityPassed',Dispatch_Check.Passed);
Task3_Result=struct('Result_Week1',Result_Week1,'Dispatch_Data',Dispatch_Data,...
 'Dispatch_Check',Dispatch_Check,'V_dispatch_ac',V_dispatch_ac,'SolverRaw',SolverRaw);
save_week1_part(script_dir,['Task3_' model_type],Task3_Result);
disp(Result_Week1);disp(Dispatch_Check);
if ProvenOptimal
    fprintf('整数最优性已在所设容差内得到证明。\n');
else
    fprintf('仅报告已验证可行解；未宣称全局最优。整数状态：%s。\n',IntegerStatus);
end

%% 9. 全时段交流电压与二阶锥松弛对照图
fig_check=figure('Color','w','Position',[100 100 900 650]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot(1:T,min(V_dispatch_ac,[],1),'-o','LineWidth',1.5);hold on;
plot(1:T,max(V_dispatch_ac,[],1),'-s','LineWidth',1.5);
yline(0.95,'--k');yline(1.05,'--k');grid on;
xlim([1 T]);ylabel('AC voltage / p.u.');legend('Minimum','Maximum','Location','best');
title(sprintf('24 h AC verification: %d, model %s',Stu_ID,model_type));
nexttile;
gaps=abs(Dispatch_Data.L.*Dispatch_Data.U(fbus,:)-Dispatch_Data.P.^2-Dispatch_Data.Q.^2);
semilogy(1:T,max(gaps,[],1)+eps,'-o','LineWidth',1.5);grid on;
xlim([1 T]);xlabel('Hour');ylabel('Maximum cone gap / p.u.');
fprintf('WEEK1_RUN_COMPLETE\n');
%% 任务三规定图件：充放电、SoC及节点电压
plot_week1_task3(Dispatch_Data,V_dispatch_ac,script_dir);
%% 手算二阶锥松弛间隙
U_hand=1.00; L_hand=0.02; P_hand=0.10; Q_hand=0.05;
g_hand=abs(U_hand*L_hand-P_hand^2-Q_hand^2);
disp(table(g_hand));
assert(abs(g_hand-0.0075)<1e-12);

function [report,V_ac] = week1_check_dispatch(d)
% 在相同负荷、光伏和储能出力下检查运行约束，并回代交流潮流。
mpc=d.mpc;nb=size(mpc.bus,1);nl=size(mpc.branch,1);T=size(d.SoC,2);
fb=mpc.branch(:,1);tb=mpc.branch(:,2);r=mpc.branch(:,3);x=mpc.branch(:,4);
af=sparse(fb,1:nl,1,nb,nl);at=sparse(tb,1:nl,1,nb,nl);
ae=sparse(d.ess_buses,1:numel(d.ess_buses),1,nb,numel(d.ess_buses));
ap=sparse(d.pv_buses,1:numel(d.pv_buses),1,nb,numel(d.pv_buses));
a0=sparse(1,1,1,nb,1);
p_res=at*(d.P-r.*d.L)-af*d.P+a0*d.Pgrid+...
    ap*d.PpvMW/d.Sbase+ae*(d.PdisMW-d.PchMW)/d.Sbase-d.PloadMW/d.Sbase;
q_res=at*(d.Q-x.*d.L)-af*d.Q+a0*d.Qgrid-d.QloadMvar/d.Sbase;
v_res=d.U(tb,:)-d.U(fb,:)+2*(r.*d.P+x.*d.Q)-(r.^2+x.^2).*d.L;
gap=abs(d.L.*d.U(fb,:)-d.P.^2-d.Q.^2);
previous=[repmat(d.SoC0,size(d.SoC,1),1),d.SoC(:,1:end-1)];
energy_res=d.EratedMWh*(d.SoC-previous)-...
    (d.eta_ch*d.PchMW-d.PdisMW/d.eta_dis)*d.dt;
report=struct('ModelType','MISOCP','Solver',d.solver,...
    'StudentID',d.Stu_ID,'Case',d.case_name,'Objective',d.model_type,...
    'SolveSeconds',d.solve_seconds,...
    'MaxSimultaneousPower_MW',max(min(d.PchMW,d.PdisMW),[],'all'),...
    'MinSoC',min(d.SoC(:)),'MaxSoC',max(d.SoC(:)),...
    'TerminalSoCError',max(abs(d.SoC(:,end)-d.SoC0)),...
    'EnergyBalanceError_MWh',max(abs(energy_res),[],'all'),...
    'PowerLimitViolation_MW',max([0;d.PchMW(:)-d.PratedMW;d.PdisMW(:)-d.PratedMW;-d.PchMW(:);-d.PdisMW(:)]),...
    'PVLimitViolation_MW',max([0;d.PpvMW(:)-d.PpvMaxMW(:);-d.PpvMW(:)]),...
    'IntegerError',max(abs(d.u(:)-round(d.u(:)))),...
    'MaxConeGap',max(gap(:)),...
    'ActiveBalanceError_pu',max(abs(p_res),[],'all'),...
    'ReactiveBalanceError_pu',max(abs(q_res),[],'all'),...
    'VoltageEquationError_pu',max(abs(v_res),[],'all'),...
    'MinModelVoltage',min(sqrt(d.U(:))),...
    'MaxModelVoltage',max(sqrt(d.U(:))),...
    'Curtailment_MWh',sum(d.PpvMaxMW-d.PpvMW,'all')*d.dt,...
    'LineLoss_MWh',d.Sbase*sum(r.*d.L,'all')*d.dt,...
    'GridCost_Yuan',sum(d.price_MWh(:).*(d.Sbase*d.Pgrid(:)))*d.dt);
V_ac=zeros(nb,T);ac_loss=zeros(1,T);ac_success=true;
opts=mpoption('verbose',0,'out.all',0);
for t=1:T
    ac_case=mpc;
    ac_case.bus(:,3)=d.PloadMW(:,t);
    ac_case.bus(:,4)=d.QloadMvar(:,t);
    ac_case.bus(d.pv_buses,3)=ac_case.bus(d.pv_buses,3)-d.PpvMW(:,t);
    ac_case.bus(d.ess_buses,3)=ac_case.bus(d.ess_buses,3)+d.PchMW(:,t)-d.PdisMW(:,t);
    ac_result=runpf(ac_case,opts);
    ac_success=ac_success && ac_result.success==1;
    V_ac(:,t)=ac_result.bus(:,8);
    ac_loss(t)=sum(real(get_losses(ac_result)));
end
report.ACConverged=ac_success;
report.MinACVoltage=min(V_ac(:));report.MaxACVoltage=max(V_ac(:));
report.MaxACVoltageError=max(abs(V_ac-sqrt(d.U)),[],'all');
voltage_error=abs(V_ac-sqrt(d.U));
[~,error_index]=max(voltage_error(:));
[report.MaxVoltageErrorNode,report.MaxVoltageErrorHour]=ind2sub(size(V_ac),error_index);
report.ACLineLoss_MWh=sum(ac_loss)*d.dt;
report.StorageConstraintsPassed=...
    report.MaxSimultaneousPower_MW<=1e-6 && report.IntegerError<=1e-6 &&...
    report.MinSoC>=d.SoCmin-1e-6 && report.MaxSoC<=d.SoCmax+1e-6 &&...
    report.TerminalSoCError<=1e-6 && report.EnergyBalanceError_MWh<=1e-6 &&...
    report.PowerLimitViolation_MW<=1e-5;
report.LinearConstraintsPassed=report.PVLimitViolation_MW<=1e-5 &&...
    report.ActiveBalanceError_pu<=1e-7 && report.ReactiveBalanceError_pu<=1e-7 &&...
    report.VoltageEquationError_pu<=1e-7;
report.ConeRelaxationExact=report.MaxConeGap<=1e-5;
report.ACVoltageWithinLimits=all(isfinite(V_ac(:))) && ac_success &&...
    report.MinACVoltage>=0.95-1e-6 && report.MaxACVoltage<=1.05+1e-6;
report.Passed=report.StorageConstraintsPassed && report.LinearConstraintsPassed &&...
    report.ConeRelaxationExact && report.ACVoltageWithinLimits;
report.TeachingCounterexample=strcmp(d.model_type,'C') &&...
    report.StorageConstraintsPassed && report.LinearConstraintsPassed && ~report.ConeRelaxationExact;
report.TeachingRequirementMet=report.Passed || report.TeachingCounterexample;
end

function plot_week1_task3(d,vac,out_dir)
%PLOT_WEEK1_TASK3 Plot dispatch and MATPOWER voltage arrays.
assert(isequal(size(vac),[size(d.mpc.bus,1),24]));
f=figure('Color','w','Position',[80 60 1100 960]);
try,theme(f,'light');catch,end
tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
labels=compose('节点%d',d.ess_buses);
nexttile;plot(1:24,(d.PdisMW-d.PchMW)','-o','LineWidth',1.4,'MarkerSize',3);
yline(0,'k:');grid on;xlim([1 24]);ylim([-0.6 0.3]);ylabel('净放电功率 / MW','FontName','SimSun');
title('充放电功率（正值放电，负值充电）','FontName','SimSun');
legend(labels,'FontName','SimSun','Location','eastoutside');
nexttile;plot(0:24,100*[repmat(d.SoC0,numel(d.ess_buses),1),d.SoC]','-o','LineWidth',1.4,'MarkerSize',3);
yline(100*d.SoCmin,'k--');yline(100*d.SoCmax,'k--');grid on;xlim([0 24]);ylim([0 100]);
ylabel('SoC / %','FontName','SimSun');title('各储能节点荷电状态','FontName','SimSun');
legend(labels,'FontName','SimSun','Location','eastoutside');
nexttile;plot(1:24,vac','Color',[0.72 0.78 0.85],'LineWidth',0.6);hold on;
[~,weak]=min(min(vac,[],2));
h1=plot(1:24,vac(weak,:),'-o','Color',[0 .35 .75],'LineWidth',1.6,'MarkerSize',3);
h2=plot(1:24,max(vac,[],1),'Color',[.9 .4 0],'LineWidth',1.6);
yline(.95,'k--');yline(1.05,'k--');grid on;xlim([1 24]);ylim([.94 1.065]);
xlabel('时段 / h','FontName','SimSun');ylabel('交流电压 / p.u.','FontName','SimSun');
title('MATPOWER回代：全部49节点电压','FontName','SimSun');
legend([h1,h2],{sprintf('最低电压所在节点%d',d.mpc.bus(weak,1)),'逐时最高电压'},'FontName','SimSun','Location','eastoutside');
sgtitle(sprintf('任务三　雷宇杰 206231313　目标%s',d.model_type),'FontName','SimSun');
set(findall(f,'Type','axes'),'FontName','Times New Roman','FontSize',14);
set(findall(f,'Type','legend'),'FontSize',13);
if ~strcmp(d.model_type,'A'),return;end
base=fullfile(out_dir,'Fig_Week1_Task3_206231313');
exportgraphics(f,[base '.png'],'Resolution',600,'BackgroundColor','white');
exportgraphics(f,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
savefig(f,[base '.fig']);
end

function save_week1_part(folder,field,data)
file=fullfile(folder,'result_summary.mat');
if isfile(file),saved=load(file);else,saved=struct;end
saved.(field)=data;
if strcmp(field,'Task3_A')
 r=data.Result_Week1;
 saved.Result=struct('StuID',r.StuID,'Stage','Week1','Obj_Value',r.ObjectiveValue,...
 'Model_Type','A','Solver',r.Solver,'SOCP_Gap',r.MaxConeGap,'Max_V_Err',r.Max_V_Err,'AC_Passed',r.FeasibilityPassed);
end
save(file,'-struct','saved','-v7');
end
