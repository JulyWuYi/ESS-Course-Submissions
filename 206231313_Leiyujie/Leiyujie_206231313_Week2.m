clear; clc; close all;
% 储能231 雷宇杰 206231313。课程模板的个人版本。
script_dir=fileparts(mfilename('fullpath'));
if isempty(script_dir),script_dir=pwd;end
project_root=script_dir;
while true
 if isfile(fullfile(project_root,'setup_env.m'))
  addpath(project_root);setup_env;break;
 end
 course_cases=fullfile(project_root,'ESS-Course-Materials','cases');
 if isfolder(course_cases),addpath(course_cases);break;end
 parent_dir=fileparts(project_root);
 if strcmp(parent_dir,project_root),break;end
 project_root=parent_dir;
end
assert(exist('runpf','file')~=0 && exist('sdpvar','file')~=0 && exist('case_nongcun','file')~=0,...
 '请先按课程说明配置MATPOWER、YALMIP、Optimization Toolbox及课程cases。');
cd(script_dir);yalmip('clear');


assert(exist('case_nongcun','file')~=0,'请先配置课程cases目录。');

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
    cand_buses = [4, 9, 14];
elseif ismember(X_digit, [3, 4, 5])
    cand_buses = [5, 11, 16];
else
    cand_buses = [6, 13, 17];
end
N_cand = length(cand_buses);

P_opt_MW  = [1, 2, 3, 4, 5]; % 离散功率档位 (MW)
E_opt_MWh = [1, 2, 3, 4, 5]; % 离散容量档位 (MWh)
K = length(E_opt_MWh);

eta_ch  = 0.90;
eta_dis = 0.95;
SoC_min = 0.20;
SoC_max = 0.80;
SoC_0   = 0.50;

c_P    = 0.08 * 1e6;   % 单位功率投资成本: 80,000 元/MW
c_E    = 0.55 * 1e6;   % 单位能量投资成本: 550,000 元/MWh
c_om   = 200 * 1e3;    % 单位功率年运维成本: 200,000 元/(MW*年)
Y_life = 20;           % 规划年限: 20 年
d_rate = 0.05;         % 贴现率: 5%
CRF    = (d_rate * (1 + d_rate)^Y_life) / ((1 + d_rate)^Y_life - 1);
D_days = 365;

V_min_Lin = 0.960;
V_max_Lin = 1.050;

Sbase = mpc.baseMVA;
N_bus = size(mpc.bus, 1);
N_br  = size(mpc.branch, 1);
fbus  = mpc.branch(:, 1);
tbus  = mpc.branch(:, 2);
r_br  = mpc.branch(:, 3);
x_br  = mpc.branch(:, 4);

Pd_pu = (mpc.bus(:, 3) * alpha / Sbase) * load_curve_pu';
Qd_pu = (mpc.bus(:, 4) * alpha / Sbase) * load_curve_pu';

z        = binvar(N_cand, 1, 'full');    % 候选节点选址二进制变量 z_i
w        = binvar(N_cand, K, 'full');    % 容量档位选择二进制变量 w_{i,k}
SoC      = sdpvar(N_cand, T, 'full');    % 储能荷电状态 SoC_{i,t}
delta    = sdpvar(N_cand, K, T, 'full'); % Big-M 辅助变量 \delta_{i,k,t} = SoC_{i,t} * w_{i,k}
P_ch_MW  = sdpvar(N_cand, T, 'full');    % 充电功率 (MW)
P_dis_MW = sdpvar(N_cand, T, 'full');    % 放电功率 (MW)
u_ch     = binvar(N_cand, T, 'full');    % 充放电互斥二进制变量

U_lin     = sdpvar(N_bus, T, 'full');    % LinDistFlow 节点电压平方 (p.u.^2)
P_br      = sdpvar(N_br,  T, 'full');    % 支路有功潮流 (p.u.)
Q_br      = sdpvar(N_br,  T, 'full');    % 支路无功潮流 (p.u.)
P_grid_pu = sdpvar(1,     T, 'full');    % 关口购电有功功率 (p.u.)
Q_grid_pu = sdpvar(1,     T, 'full');    % 关口注入无功功率 (p.u.)

Constraints = [];
P_rated     = w * P_opt_MW';
E_rated     = w * E_opt_MWh';
P_rated_MW  = P_rated;
E_rated_MWh = E_rated;
P_max_bound = max(P_opt_MW);

for idx = 1:N_cand
    Constraints = [Constraints, sum(w(idx, :)) == z(idx)];
    Constraints = [Constraints, SoC_min <= SoC(idx, :) <= SoC_max];
    Constraints = [Constraints, SoC(idx, T) == SoC_0];
    
    Constraints = [Constraints, 0 <= P_ch_MW(idx, :)  <= P_rated_MW(idx)];
    Constraints = [Constraints, 0 <= P_dis_MW(idx, :) <= P_rated_MW(idx)];
    Constraints = [Constraints, P_ch_MW(idx, :)  <= u_ch(idx, :) * P_max_bound];
    Constraints = [Constraints, P_dis_MW(idx, :) <= (1 - u_ch(idx, :)) * P_max_bound];
    
    % 【完成 1】双线性乘积项：双线性乘积项 \delta_{i,k,t} = SoC_{i,t} * w_{i,k} 的四条大M法线性化不等式：
    for k = 1:K
        d_vec = reshape(delta(idx, k, :), 1, T);
        Constraints = [Constraints, SoC_min * w(idx, k) <= d_vec <= SoC_max * w(idx, k)];
        Constraints = [Constraints, SoC(idx, :) - SoC_max * (1 - w(idx, k)) <= d_vec <= ...
            SoC(idx, :) - SoC_min * (1 - w(idx, k))];
    end
    
    % 【完成 2】线性化时序能量：线性化后的时序能量守恒递推方程：
    for t = 1:T
        if t == 1
            E_prev = E_opt_MWh * (SoC_0 * w(idx, :))';
        else
            E_prev = E_opt_MWh * reshape(delta(idx, :, t - 1), K, 1);
        end
        E_curr = E_opt_MWh * reshape(delta(idx, :, t), K, 1);
        Constraints = [Constraints, E_curr - E_prev == ...
            eta_ch * P_ch_MW(idx, t) - P_dis_MW(idx, t) / eta_dis];
    end
end

A_from  = sparse(fbus, 1:N_br, 1, N_bus, N_br);
A_to    = sparse(tbus, 1:N_br, 1, N_bus, N_br);
A_cand  = sparse(cand_buses, 1:N_cand, 1, N_bus, N_cand);
A_slack = sparse(1, 1, 1, N_bus, 1);

Constraints = [Constraints, A_to * P_br - A_from * P_br + A_slack * P_grid_pu + ...
    A_cand * ((P_dis_MW - P_ch_MW) / Sbase) == Pd_pu];
Constraints = [Constraints, A_to * Q_br - A_from * Q_br + A_slack * Q_grid_pu == Qd_pu];
Constraints = [Constraints, U_lin(1, :) == 1.0];
Constraints = [Constraints, (V_min_Lin^2) <= U_lin(2:end, :) <= (V_max_Lin^2)];

% 【完成 3】LinDistFlow： LinDistFlow 支路电压平方降落方程：
for b = 1:N_br
    Constraints = [Constraints, U_lin(tbus(b), :) == U_lin(fbus(b), :) - ...
        2 * (r_br(b) * P_br(b, :) + x_br(b) * Q_br(b, :))];
end

C_inv_ann  = CRF * sum(c_P * P_rated_MW + c_E * E_rated_MWh);
C_om_ann   = sum(c_om * P_rated_MW);
C_grid_ann = D_days * sum((Sbase * P_grid_pu') .* price_MWh);
Objective  = C_inv_ann + C_om_ann + C_grid_ann;

native_optim_dir = fullfile(matlabroot, 'toolbox', 'optim', 'optim');
if exist(native_optim_dir, 'dir')
    addpath(native_optim_dir, '-begin');
end
solver_name = 'intlinprog';
ops = sdpsettings('solver', solver_name, 'verbose', 1, 'savesolveroutput', 1);
ops.intlinprog.RelativeGapTolerance = 1e-6;
ops.intlinprog.AbsoluteGapTolerance = 1e-8;
solve_clock=tic;
sol = optimize(Constraints, Objective, ops);
solve_seconds=toc(solve_clock);

if sol.problem ~= 0
    error('MILP 选址定容求解失败: %s', yalmiperror(sol.problem));
end

V_lin_pu      = sqrt(value(U_lin));
V_ac_base     = zeros(N_bus, T);
V_ac_ess      = zeros(N_bus, T);
Ploss_ac_base = zeros(T, 1);
Ploss_ac_ess  = zeros(T, 1);
P_net_ess_MW  = value(P_dis_MW - P_ch_MW);
mpopt = mpoption('verbose', 0, 'out.all', 0);

for t = 1:T
    mpc_b = mpc;
    mpc_b.bus(:, 3) = Pd_pu(:, t) * Sbase;
    mpc_b.bus(:, 4) = Qd_pu(:, t) * Sbase;
    res_b = runpf(mpc_b, mpopt);
    assert(res_b.success==1,'基态第%d时段交流潮流未收敛。',t);
    V_ac_base(:, t)  = res_b.bus(:, 8);
    Ploss_ac_base(t) = sum(real(get_losses(res_b)));
    
    mpc_e = mpc_b;
    mpc_e.bus(cand_buses, 3) = mpc_e.bus(cand_buses, 3) - P_net_ess_MW(:, t);
    res_e = runpf(mpc_e, mpopt);
    assert(res_e.success==1,'储能第%d时段交流潮流未收敛。',t);
    V_ac_ess(:, t)  = res_e.bus(:, 8);
    Ploss_ac_ess(t) = sum(real(get_losses(res_e)));
end

B_loss_ann           = D_days * sum((Ploss_ac_base - Ploss_ac_ess) .* price_MWh);
max_voltage_error_pu = max(max(abs(V_lin_pu - V_ac_ess)));

fprintf('\n================ 第二周：储能选址定容与交流潮流验证报告 ================\n');
fprintf('学号: %d | 算例: %s | 求解器: %s\n', Stu_ID, case_name, solver_name);
for idx = 1:N_cand
    fprintf('候选节点 %2d: 选址状态 z=%d | 额定功率 = %d MW | 额定容量 = %d MWh\n', ...
        cand_buses(idx), round(value(z(idx))), round(value(P_rated_MW(idx))), round(value(E_rated_MWh(idx))));
end
fprintf('年化投资成本: %.2f 万元 | 年运维成本: %.2f 万元 | 事后核算年降损收益: %.2f 万元\n', ...
    value(C_inv_ann) / 1e4, value(C_om_ann) / 1e4, B_loss_ann / 1e4);
fprintf('LinDistFlow 最低电压: %.4f p.u. | 交流潮流真实最低电压: %.4f p.u. | 最大电压偏差: %.4f p.u.\n', ...
    min(V_lin_pu(:)), min(V_ac_ess(:)), max_voltage_error_pu);
fprintf('========================================================================\n');

Result.StuID       = Stu_ID;                   % 学生学号
Result.Obj_Value   = value(Objective);         % 年综合费用最优值 (元)
Result.Opt_Sites   = find(round(value(z))==1); % 选址候选节点集合中的序号 (1..3，对应母线 cand_buses(Result.Opt_Sites))
Result.Opt_P_rated = value(P_rated);           % 各候选节点配置功率 (MW)
Result.Opt_E_rated = value(E_rated);           % 各候选节点配置能量 (MWh)
Result.Max_V_Err   = max_voltage_error_pu;     % LinDistFlow 与 AC 潮流最大电压误差 (p.u.)
save_week2_summary(script_dir,Result);
fprintf('[OK] 标准化评测结果已自动保存至当前目录的 result_summary.mat\n');

[~, worst_bus] = min(min(V_ac_base, [], 2));
fig1 = figure('Name', '薄弱节点24h电压幅值对比与线性化误差校验', 'Color', 'w', 'Position', [150, 150, 740, 440]);
if isprop(fig1, 'Theme'), theme(fig1, 'light'); end
p1 = plot(1:T, V_ac_base(worst_bus, :), '--k', 'LineWidth', 1.5); hold on;
p2 = plot(1:T, V_lin_pu(worst_bus, :), '-b^', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', [0.00, 0.45, 0.74]);
p3 = plot(1:T, V_ac_ess(worst_bus, :), '-ro', 'LineWidth', 1.7, 'MarkerSize', 5, 'MarkerFaceColor', [0.85, 0.20, 0.15], 'Color', [0.85, 0.20, 0.15]);
p4 = yline(0.95, ':r', 'LineWidth', 1.4);

y_min_plot = min([V_ac_base(worst_bus, :), V_ac_ess(worst_bus, :), 0.948]) - 0.004;
y_max_plot = max([V_lin_pu(worst_bus, :), V_ac_ess(worst_bus, :), 0.980]) + 0.022;
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontName', 'Times New Roman', 'FontSize', 10.5, 'LineWidth', 1.0, ...
    'Box', 'on', 'TickDir', 'in', 'XTick', 0:4:24, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'GridColor', [0.4, 0.4, 0.4]);
xlim([1, T]);
ylim([y_min_plot, y_max_plot]);
xlabel('\fontname{宋体}时间 \fontname{Times New Roman}{\itt}\rm{/h}', 'FontSize', 12, 'Color', 'k', 'Interpreter', 'tex');
ylabel('\fontname{宋体}节点电压幅值 \fontname{Times New Roman}{\itV}\rm{/(p.u.)}', 'FontSize', 12, 'Color', 'k', 'Interpreter', 'tex');
title(sprintf('\\fontname{宋体}薄弱节点 \\fontname{Times New Roman}%d \\fontname{宋体}的 \\fontname{Times New Roman}24 h \\fontname{宋体}电压对比 (\\fontname{Times New Roman}\\Delta{\\itV}_{\\rm{max}} = %.4f p.u.\\fontname{宋体})', ...
    worst_bus, max_voltage_error_pu), 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'k', 'Interpreter', 'tex');
leg1 = legend([p1, p2, p3, p4], { ...
    '\fontname{宋体}未配置储能 (\fontname{Times New Roman}AC \fontname{宋体}潮流 \fontname{Times New Roman}{\itV}_{\rm{AC}}^{\rm{base}}\fontname{宋体})', ...
    '\fontname{宋体}配置储能后 (\fontname{Times New Roman}LinDistFlow \fontname{Times New Roman}{\itV}_{\rm{Lin}}^{\rm{ess}}\fontname{宋体})', ...
    '\fontname{宋体}配置储能后 (\fontname{Times New Roman}AC \fontname{宋体}潮流 \fontname{Times New Roman}{\itV}_{\rm{AC}}^{\rm{ess}}\fontname{宋体})', ...
    '\fontname{宋体}电压安全下限 \fontname{Times New Roman}{\itV}_{\rm{min}} = 0.95 p.u.'}, ...
    'Location', 'north', 'NumColumns', 2, 'FontSize', 10, 'Interpreter', 'tex', 'Box', 'on');
set(leg1, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
grid on;

% 默认自动导出高质量白底 PNG (600 dpi)、矢量 PDF 与 MATLAB 源文件 FIG
exportgraphics(fig1, sprintf('Fig_Week2_%d.png', Stu_ID), 'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(fig1, sprintf('Fig_Week2_%d.pdf', Stu_ID), 'ContentType', 'vector', 'BackgroundColor', 'white');
savefig(fig1, sprintf('Fig_Week2_%d.fig', Stu_ID));

%% 独立数值校核与逐时结果归档
Plan=struct('StuID',Stu_ID,'Name','雷宇杰','Class','储能231',...
 'case_name',case_name,'alpha',alpha,'mpc',mpc,'Sbase',Sbase,'T',T,'dt',1,...
 'cand_buses',cand_buses,'P_options',P_opt_MW,'E_options',E_opt_MWh,...
 'z',value(z),'w',value(w),'SoC',value(SoC),'delta',value(delta),...
 'Pch_MW',value(P_ch_MW),'Pdis_MW',value(P_dis_MW),'u',value(u_ch),...
 'P_rated_MW',value(P_rated),'E_rated_MWh',value(E_rated),...
 'eta_ch',eta_ch,'eta_dis',eta_dis,'SoC_min',SoC_min,'SoC_max',SoC_max,'SoC0',SoC_0,...
 'U_lin',value(U_lin),'P_br',value(P_br),'Q_br',value(Q_br),...
 'Pgrid_pu',value(P_grid_pu),'Qgrid_pu',value(Q_grid_pu),'Pd_pu',Pd_pu,'Qd_pu',Qd_pu,...
 'V_min_Lin',V_min_Lin,'V_max_Lin',V_max_Lin,'price',price_MWh(:),...
 'CRF',CRF,'d_rate',d_rate,'Y_life',Y_life,'D_days',D_days,...
 'c_P',c_P,'c_E',c_E,'c_om',c_om,'C_inv_ann',value(C_inv_ann),...
 'C_om_ann',value(C_om_ann),'C_grid_ann',value(C_grid_ann),'Objective',value(Objective),...
 'V_ac_base',V_ac_base,'V_ac_ess',V_ac_ess,'Ploss_ac_base',Ploss_ac_base,...
 'Ploss_ac_ess',Ploss_ac_ess,'B_loss_ann',B_loss_ann,...
 'solver_name',solver_name,'solve_seconds',solve_seconds,'solver_output',sol.solveroutput,'YALMIP_problem',sol.problem);
Check=check_week2_plan(Plan);
assert(Check.AllConstraintsPassed,'规划解的独立数值校验未通过。');
assert(Check.ACVoltagePassed,'储能方案的真实交流电压越限。');
audit_dir=week2_audit_directory(project_root);
save(fullfile(audit_dir,'Week2_Full_Result.mat'),'Plan','Check','Result','-v7');
save_week2_summary(script_dir,Result,Plan,Check);
facts=week2_facts(Plan,Check,Result);
write_week2_json(fullfile(audit_dir,'Week2_Facts.json'),facts);
hours=(1:T)';
schedule=table(hours,price_MWh(:),Sbase*value(P_grid_pu(:)),min(V_ac_base,[],1)',...
 min(V_ac_ess,[],1)',max(V_ac_ess,[],1)',min(V_lin_pu,[],1)',...
 Ploss_ac_base(:),Ploss_ac_ess(:),'VariableNames',...
 {'Hour','Price_YuanPerMWh','LinGrid_MW','MinACBase_pu','MinACStorage_pu','MaxACStorage_pu','MinLin_pu','ACBaseLoss_MW','ACStorageLoss_MW'});
for site=1:N_cand
 schedule.(sprintf('Node%d_Charge_MW',cand_buses(site)))=value(P_ch_MW(site,:))';
 schedule.(sprintf('Node%d_Discharge_MW',cand_buses(site)))=value(P_dis_MW(site,:))';
 schedule.(sprintf('Node%d_SoC',cand_buses(site)))=value(SoC(site,:))';
end
writetable(schedule,fullfile(audit_dir,'Week2_Timeseries.csv'));
disp(facts);fprintf('WEEK2_NUMERICAL_AUDIT_PASSED\n');

function save_week2_summary(folder,Result,Plan,Check)
file=fullfile(folder,'result_summary.mat');
if isfile(file),saved=load(file);else,saved=struct;end
saved.Result=Result;
if nargin>2,saved.Week2_Result=struct('Plan',Plan,'Check',Check,'Result',Result);end
save(file,'-struct','saved','-v7');
end

function folder=week2_audit_directory(project_root)
folder=getenv('ESS_WEEK2_AUDIT_DIR');
if isempty(folder)
 if isfolder(fullfile(project_root,'logs'))
  folder=fullfile(project_root,'logs','脚本审计');
 else
  folder=fullfile(tempdir,'ESS_Week2_206231313');
 end
end
if ~isfolder(folder),mkdir(folder);end
end

function c=check_week2_plan(d)
nb=size(d.mpc.bus,1);nl=size(d.mpc.branch,1);nc=numel(d.cand_buses);T=d.T;
fb=d.mpc.branch(:,1);tb=d.mpc.branch(:,2);r=d.mpc.branch(:,3);x=d.mpc.branch(:,4);
af=sparse(fb,1:nl,1,nb,nl);at=sparse(tb,1:nl,1,nb,nl);
ae=sparse(d.cand_buses,1:nc,1,nb,nc);a0=sparse(1,1,1,nb,1);
pr=at*d.P_br-af*d.P_br+a0*d.Pgrid_pu+ae*(d.Pdis_MW-d.Pch_MW)/d.Sbase-d.Pd_pu;
qr=at*d.Q_br-af*d.Q_br+a0*d.Qgrid_pu-d.Qd_pu;
vr=d.U_lin(tb,:)-d.U_lin(fb,:)+2*(r.*d.P_br+x.*d.Q_br);
product=zeros(size(d.delta));viol=0;
for i=1:nc
 for k=1:numel(d.E_options)
  dv=reshape(d.delta(i,k,:),1,T);wv=d.w(i,k);
  product(i,k,:)=reshape(d.SoC(i,:)*wv,1,1,T);
  viol=max([viol,max(d.SoC_min*wv-dv),max(dv-d.SoC_max*wv),...
   max(d.SoC(i,:)-d.SoC_max*(1-wv)-dv),max(dv-d.SoC(i,:)+d.SoC_min*(1-wv))]);
 end
end
prev=[repmat(d.SoC0,nc,1),d.SoC(:,1:end-1)];
energy=d.E_rated_MWh.*(d.SoC-prev)-(d.eta_ch*d.Pch_MW-d.Pdis_MW/d.eta_dis)*d.dt;
limit=max([0;d.Pch_MW(:)-repelem(d.P_rated_MW,T);d.Pdis_MW(:)-repelem(d.P_rated_MW,T)]);
% Arrays are column-major, so use a full repeated rating matrix.
ratings=repmat(d.P_rated_MW,1,T);
limit=max([0;reshape(d.Pch_MW-ratings,[],1);reshape(d.Pdis_MW-ratings,[],1);-d.Pch_MW(:);-d.Pdis_MW(:)]);
overlap=max(min(d.Pch_MW,d.Pdis_MW),[],'all');
integer=max(abs([d.z(:);d.w(:);d.u(:)]-round([d.z(:);d.w(:);d.u(:)])));
choice=max(abs(sum(d.w,2)-d.z));
rating_error=max([abs(d.P_rated_MW-d.w*d.P_options(:));abs(d.E_rated_MWh-d.w*d.E_options(:))]);
vb=zeros(nb,T);ve=vb;lb=zeros(T,1);le=lb;gridb=lb;gride=lb;converged=true;
op=mpoption('verbose',0,'out.all',0);
for t=1:T
 m=d.mpc;m.bus(:,3)=d.Pd_pu(:,t)*d.Sbase;m.bus(:,4)=d.Qd_pu(:,t)*d.Sbase;
 a=runpf(m,op);converged=converged && a.success==1;vb(:,t)=a.bus(:,8);lb(t)=sum(real(get_losses(a)));gridb(t)=sum(a.gen(:,2));
 m.bus(d.cand_buses,3)=m.bus(d.cand_buses,3)+d.Pch_MW(:,t)-d.Pdis_MW(:,t);
 a=runpf(m,op);converged=converged && a.success==1;ve(:,t)=a.bus(:,8);le(t)=sum(real(get_losses(a)));gride(t)=sum(a.gen(:,2));
end
ac_array=max([abs(vb(:)-d.V_ac_base(:));abs(ve(:)-d.V_ac_ess(:));abs(lb-d.Ploss_ac_base);abs(le-d.Ploss_ac_ess)]);
dv=abs(sqrt(d.U_lin)-ve);[err,idx]=max(dv(:));[node,hour]=ind2sub(size(dv),idx);
[minbase,idx]=min(vb(:));[base_node,base_hour]=ind2sub(size(vb),idx);
[miness,idx]=min(ve(:));[ess_node,ess_hour]=ind2sub(size(ve),idx);
c=struct('ActiveBalanceError_pu',max(abs(pr),[],'all'),'ReactiveBalanceError_pu',max(abs(qr),[],'all'),...
 'LinVoltageEquationError_pu',max(abs(vr),[],'all'),'BigMProductError',max(abs(d.delta-product),[],'all'),...
 'BigMInequalityViolation',max(viol,0),'EnergyBalanceError_MWh',max(abs(energy),[],'all'),...
 'PowerLimitViolation_MW',limit,'MaxSimultaneousPower_MW',overlap,'IntegerError',integer,...
 'ChoiceEqualityError',choice,'RatingEqualityError',rating_error,...
 'TerminalSoCError',max(abs(d.SoC(:,end)-d.SoC0)),'MinSoC',min(d.SoC(:)),'MaxSoC',max(d.SoC(:)),...
 'ACConverged',converged,'IndependentACArrayError',ac_array,'MaxVoltageError',err,...
 'MaxVoltageErrorNode',d.mpc.bus(node,1),'MaxVoltageErrorHour',hour,...
 'MinBaseACVoltage',minbase,'MinBaseACNode',d.mpc.bus(base_node,1),'MinBaseACHour',base_hour,...
 'MinStorageACVoltage',miness,'MinStorageACNode',d.mpc.bus(ess_node,1),'MinStorageACHour',ess_hour,...
 'MaxStorageACVoltage',max(ve(:)),'MinLinVoltage',min(sqrt(d.U_lin),[],'all'),...
 'BaseDailyLoss_MWh',sum(lb)*d.dt,'StorageDailyLoss_MWh',sum(le)*d.dt,...
 'BaseACGridAnnualCost_Yuan',d.D_days*sum(gridb.*d.price)*d.dt,...
 'StorageACGridAnnualCost_Yuan',d.D_days*sum(gride.*d.price)*d.dt,...
 'ACLossBenefitAnnual_Yuan',d.D_days*sum((lb-le).*d.price)*d.dt);
c.StorageACAnnualTotalCost_Yuan=d.C_inv_ann+d.C_om_ann+c.StorageACGridAnnualCost_Yuan;
c.ACNetAnnualBenefit_Yuan=c.BaseACGridAnnualCost_Yuan-c.StorageACAnnualTotalCost_Yuan;
c.ACVoltagePassed=converged && miness>=.95-1e-6 && max(ve(:))<=1.05+1e-6;
c.AllConstraintsPassed=max([c.ActiveBalanceError_pu,c.ReactiveBalanceError_pu,c.LinVoltageEquationError_pu])<1e-7 && ...
 max([c.BigMProductError,c.BigMInequalityViolation,c.EnergyBalanceError_MWh,c.IntegerError,c.ChoiceEqualityError,c.RatingEqualityError,c.TerminalSoCError])<1e-5 && ...
 c.PowerLimitViolation_MW<1e-6 && overlap<1e-6 && c.MinSoC>=d.SoC_min-1e-6 && c.MaxSoC<=d.SoC_max+1e-6 && ac_array<1e-9;
assert(abs(c.ACLossBenefitAnnual_Yuan-d.B_loss_ann)<1e-6);
end

function f=week2_facts(d,c,Result)
f=struct('student_id',d.StuID,'name',d.Name,'class',d.Class,'case_name',d.case_name,'alpha',d.alpha,...
 'solver',d.solver_name,'solve_seconds',d.solve_seconds,'YALMIP_problem',d.YALMIP_problem,...
 'candidate_buses',d.cand_buses,'selected_indices',Result.Opt_Sites,'selected_buses',d.cand_buses(Result.Opt_Sites),...
 'z',d.z,'w',d.w,'power_MW',d.P_rated_MW,'energy_MWh',d.E_rated_MWh,...
 'CRF',d.CRF,'discount_rate',d.d_rate,'lifetime_years',d.Y_life,'days',d.D_days,...
 'c_P',d.c_P,'c_E',d.c_E,'c_om',d.c_om,'V_min_Lin',d.V_min_Lin,...
 'annual_investment_yuan',d.C_inv_ann,'annual_om_yuan',d.C_om_ann,...
 'annual_lin_grid_yuan',d.C_grid_ann,'annual_objective_yuan',d.Objective,...
 'daily_charge_MWh',sum(d.Pch_MW,2)*d.dt,'daily_discharge_MWh',sum(d.Pdis_MW,2)*d.dt,...
 'check',c,'proven_optimal_within_solver_tolerance',d.YALMIP_problem==0);
f.solver_certificate=rmfield(d.solver_output,'x');
end

function write_week2_json(file,obj)
f=fopen(file,'w','n','UTF-8');fprintf(f,'%s',jsonencode(obj,PrettyPrint=true));fclose(f);
end
