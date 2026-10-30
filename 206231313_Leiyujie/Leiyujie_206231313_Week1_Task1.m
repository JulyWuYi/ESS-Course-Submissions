%% 任务一：24小时交流潮流及增载对照
% 储能231 雷宇杰 206231313。节点4/9/14各0.8 MW光伏，无储能。 运行前按课程要求配置依赖；计算与绘图函数均在本文件末尾。

task_dir=fileparts(mfilename('fullpath'));
if isempty(task_dir),task_dir=pwd;end
project_root=task_dir;
while true
 if isfile(fullfile(project_root,'setup_env.m'))
  addpath(project_root);setup_env;break;
 end
 candidate=fullfile(project_root,'ESS-Course-Materials','cases');
 if isfolder(candidate),addpath(candidate);break;end
 parent_dir=fileparts(project_root);
 if strcmp(parent_dir,project_root),break;end
 project_root=parent_dir;
end
assert(exist('runpf','file')~=0 && exist('case_nongcun','file')~=0 && exist('sdpvar','file')~=0,...
 '请先按课程说明配置MATPOWER、YALMIP、Optimization Toolbox和课程cases。');
close all;

Task1_Result=task1_experiment(task_dir);
disp(Task1_Result.summary);
% 电压越限按节点-时段点统计；时段编号1—24，每段1小时。
fprintf('TASK1_COMPLETE: 48 converged power flows.\n');
function result=task1_experiment(out_dir)
%TASK1_EXPERIMENT Recompute personal baseline/stress AC cases from teacher data.
student_id=206231313;alpha=1+0.002*mod(student_id,100);dt=1;
mpc=loadcase('case_nongcun');
profiles=load('profiles_24h.mat');
assert(size(mpc.bus,1)==49 && size(mpc.branch,1)==48 && profiles.T==24);
assert(abs(alpha-1.026)<1e-12 && all(profiles.pv_curve_pu>=0));
assert(isequal(mpc.bus(:,1),(1:49)'),'This experiment expects consecutive bus IDs.');
base=solve_case(mpc,profiles,student_id,alpha,1.00,dt);
stress=solve_case(mpc,profiles,student_id,alpha,1.10,dt);
verify_task1(base,stress);

summary=table(["Baseline";"Load_x1p10"],[1;1.10],[base.alpha;stress.alpha], ...
    [base.min_voltage;stress.min_voltage],[base.min_bus;stress.min_bus], ...
    [base.min_hour;stress.min_hour],[base.max_voltage;stress.max_voltage], ...
    [base.max_bus;stress.max_bus],[base.max_hour;stress.max_hour], ...
    [base.daily_loss_mwh;stress.daily_loss_mwh], ...
    [base.low_voltage_count;stress.low_voltage_count], ...
    [base.high_voltage_count;stress.high_voltage_count], ...
    'VariableNames',{'Scenario','AdditionalLoadMultiplier','EffectiveAlpha', ...
    'MinVoltage_pu','MinBus','MinHour','MaxVoltage_pu','MaxBus','MaxHour', ...
    'DailyLoss_MWh','LowVoltagePointCount','HighVoltagePointCount'});
result=struct('student_id',student_id,'name','雷宇杰','class','储能231', ...
    'source_case','case_nongcun.m','matlab_version',version, ...
    'matpower_version',mpver,'generated_at',char(datetime('now')), ...
    'base',base,'stress',stress,'summary',summary);
metrics=struct('student_id',student_id,'alpha',alpha,'min_voltage',base.min_voltage, ...
    'min_bus',base.min_bus,'min_hour',base.min_hour,'max_voltage',base.max_voltage, ...
    'daily_loss_mwh',base.daily_loss_mwh,'low_voltage_count',base.low_voltage_count, ...
    'high_voltage_count',base.high_voltage_count,'stress_min_voltage',stress.min_voltage, ...
    'stress_daily_loss_mwh',stress.daily_loss_mwh,'stress_min_bus',stress.min_bus, ...
    'stress_min_hour',stress.min_hour,'stress_alpha',stress.alpha, ...
    'stress_max_voltage',stress.max_voltage, ...
    'stress_low_voltage_count',stress.low_voltage_count, ...
    'stress_high_voltage_count',stress.high_voltage_count, ...
    'all_48_pf_converged',all(base.success==1)&&all(stress.success==1), ...
    'max_active_balance_error_mw',max(abs([base.active_balance_error_mw;stress.active_balance_error_mw])), ...
    'voltage_limit_count_definition','49 buses times 24 hours; strict limits 0.95 and 1.05 pu', ...
    'matlab_version',version,'generated_at',result.generated_at);
result.metrics=metrics;
save_week1_part(out_dir,'Task1_Result',result);
plot_task1_voltage(base,stress,out_dir);
fprintf('TASK1 baseline: alpha=%.6f minV=%.9f bus=%d hour=%d maxV=%.9f loss=%.9f MWh lower=%d upper=%d\n', ...
    base.alpha,base.min_voltage,base.min_bus,base.min_hour,base.max_voltage, ...
    base.daily_loss_mwh,base.low_voltage_count,base.high_voltage_count);
fprintf('TASK1 loadx1.10: alpha=%.6f minV=%.9f bus=%d hour=%d maxV=%.9f loss=%.9f MWh lower=%d upper=%d\n', ...
    stress.alpha,stress.min_voltage,stress.min_bus,stress.min_hour,stress.max_voltage, ...
    stress.daily_loss_mwh,stress.low_voltage_count,stress.high_voltage_count);
fprintf('TASK1 power-balance max residual=%.3e MW\n',metrics.max_active_balance_error_mw);
end

function r=solve_case(mpc,profiles,student_id,alpha,multiplier,dt)
T=profiles.T;N=size(mpc.bus,1);pv_buses=[4 9 14];
effective_alpha=alpha*multiplier;
V=zeros(N,T);loss=zeros(T,1);success=zeros(T,1);balance=zeros(T,1);
grid_mw=zeros(T,1);net_load=zeros(T,1);pv_mw=zeros(T,1);
options=mpoption('verbose',0,'out.all',0);
for t=1:T
    hourly=mpc;
    hourly.bus(:,3)=effective_alpha*mpc.bus(:,3)*profiles.load_curve_pu(t);
    hourly.bus(:,4)=effective_alpha*mpc.bus(:,4)*profiles.load_curve_pu(t);
    hourly.bus(pv_buses,3)=hourly.bus(pv_buses,3)-0.8*profiles.pv_curve_pu(t);
    ac=runpf(hourly,options);
    assert(ac.success==1,'Task1:PF','Load multiplier %.2f, hour %d did not converge.',multiplier,t);
    V(:,t)=ac.bus(:,8);loss(t)=sum(real(get_losses(ac)));success(t)=ac.success;
    grid_mw(t)=sum(ac.gen(ac.gen(:,8)>0,2));net_load(t)=sum(hourly.bus(:,3));
    pv_mw(t)=numel(pv_buses)*0.8*profiles.pv_curve_pu(t);
    balance(t)=grid_mw(t)-net_load(t)-loss(t);
    assert(abs(balance(t))<1e-6,'Task1:HourlyBalance','Hour %d active-power balance failed.',t);
end
[minimum,kmin]=min(V(:));[row_min,hmin]=ind2sub(size(V),kmin);
[maximum,kmax]=max(V(:));[row_max,hmax]=ind2sub(size(V),kmax);
r=struct('student_id',student_id,'case_name','case_nongcun','alpha',effective_alpha, ...
    'load_multiplier',multiplier,'dt',dt,'pv_buses',pv_buses,'pv_peak_per_node_mw',0.8, ...
    'bus_ids',mpc.bus(:,1),'voltage',V,'loss_mw',loss,'success',success, ...
    'grid_mw',grid_mw,'net_load_mw',net_load,'pv_mw',pv_mw, ...
    'active_balance_error_mw',balance,'min_voltage',minimum,'min_bus',mpc.bus(row_min,1), ...
    'min_hour',hmin,'max_voltage',maximum,'max_bus',mpc.bus(row_max,1),'max_hour',hmax, ...
    'daily_loss_mwh',sum(loss)*dt,'low_voltage_count',nnz(V<0.95), ...
    'high_voltage_count',nnz(V>1.05));
end

function verify_task1(r, stress)
% Numerical acceptance checks, independent of plots and saved summaries.
assert(r.student_id==206231313 && strcmp(r.case_name,'case_nongcun'), ...
    'Task1:Identity','Wrong student or feeder.');
assert(isequal(size(r.voltage),[49 24]) && all(isfinite(r.voltage),'all') ...
    && all(r.voltage>0,'all'),'Task1:Voltage','Invalid voltage matrix.');
assert(numel(r.success)==24 && all(r.success==1), ...
    'Task1:Convergence','Every AC solve must converge.');
assert(max(abs(r.active_balance_error_mw))<1e-6, ...
    'Task1:Balance','AC generation/load/loss balance failed.');
assert(all(r.loss_mw>=-1e-9) && r.daily_loss_mwh>0, ...
    'Task1:Loss','Losses must be physically nonnegative.');
assert(abs(r.daily_loss_mwh-sum(r.loss_mw)*r.dt)<1e-10, ...
    'Task1:LossIntegral','MW-to-MWh integration is inconsistent.');
[vmin,ix]=min(r.voltage(:));
[row,hour]=ind2sub(size(r.voltage),ix);
assert(abs(vmin-r.min_voltage)<1e-12 && r.min_bus==r.bus_ids(row) ...
    && r.min_hour==hour,'Task1:Minimum','Minimum coordinates are inconsistent.');
assert(r.low_voltage_count==nnz(r.voltage<0.95) ...
    && r.high_voltage_count==nnz(r.voltage>1.05), ...
    'Task1:Limits','Violation counts must cover all bus-hour points.');
assert(abs(r.alpha-1.026*r.load_multiplier)<1e-12, ...
    'Task1:Alpha','Personal and stress scaling were mixed.');
if nargin>1
    verify_task1(stress);
    assert(abs(r.load_multiplier-1)<1e-12 ...
        && abs(stress.load_multiplier-1.10)<1e-12, ...
        'Task1:Scenarios','The requested 1.00/1.10 scenarios are missing.');
    assert(stress.min_voltage<r.min_voltage ...
        && stress.daily_loss_mwh>r.daily_loss_mwh, ...
        'Task1:Sensitivity','The heavier-load scenario needs investigation.');
end
end

function plot_task1_voltage(base,stress,out_dir)
% 主图沿用Demo01的三维电压曲面，增载对照在另一个窗口显示。
T=size(base.voltage,2);N_bus=numel(base.bus_ids);
V_bus_24h=base.voltage;
[min_V,linear_idx]=min(V_bus_24h(:));
[min_bus,min_hour]=ind2sub(size(V_bus_24h),linear_idx);
fig1 = figure('Name', '配电网24h基态节点电压分布', 'Color', 'w', 'Position', [100, 100, 720, 440]);
if isprop(fig1, 'Theme'), theme(fig1, 'light'); end
surface_handle=surf(1:T, 1:N_bus, V_bus_24h, 'EdgeAlpha', 0.45);
shading interp;
hold on;
plot3(min_hour, min_bus, min_V, 'rp', 'MarkerSize', 11, 'MarkerFaceColor', 'r');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
    'FontName', 'Times New Roman', 'FontSize', 10.5, 'LineWidth', 1.0, ...
    'Box', 'on', 'TickDir', 'in', 'GridLineStyle', '--', 'GridAlpha', 0.3, 'GridColor', [0.4, 0.4, 0.4]);
xlim([1, T]);
ylim([1, N_bus]);
xticks(0:4:24);
xlabel('\fontname{宋体}时间 \fontname{Times New Roman}{\itt}\rm{/h}', 'FontSize', 12, 'Color', 'k', 'Interpreter', 'tex');
ylabel('\fontname{宋体}节点编号 \fontname{Times New Roman}{\iti}', 'FontSize', 12, 'Color', 'k', 'Interpreter', 'tex');
zlabel('\fontname{宋体}节点电压幅值 \fontname{Times New Roman}{\itV}\rm{/(p.u.)}', 'FontSize', 12, 'Color', 'k', 'Interpreter', 'tex');
title(sprintf('\\fontname{宋体}配电网 24 h 节点电压分布 \\fontname{Times New Roman}({\\itV}_{\\rm{min}} = %.4f p.u. @ Bus %d, {\\itt} = %d h)', ...
    min_V, min_bus, min_hour), 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'k', 'Interpreter', 'tex');
cb = colorbar;
cb.FontName = 'Times New Roman';
cb.FontSize = 10.5;
cb.Color = 'k';
cb.Label.String = '\fontname{宋体}电压幅值 \fontname{Times New Roman}{\itV}\rm{/(p.u.)}';
cb.Label.Interpreter = 'tex';
cb.Label.FontSize = 11;
cb.Label.Color = 'k';
view(-38, 28);
grid on;

assert(isequal(surface_handle.ZData,base.voltage));
stem=fullfile(out_dir,'Fig_Week1_Task1_206231313');
exportgraphics(fig1,[stem '.png'],'Resolution',600,'BackgroundColor','white');
% 与Demo01一致，三维渐变曲面在PDF内以600 dpi保存。
exportgraphics(fig1,[stem '.pdf'],'ContentType','image','Resolution',600,'BackgroundColor','white');
savefig(fig1,[stem '.fig']);
plot_task1_comparison(base,stress);
end

function plot_task1_comparison(base,stress)
% Replot the original 49-by-24 samples without changing any numerical result.
assert(isequal(size(base.voltage),[49 24]) && isequal(size(stress.voltage),[49 24]));
assert(isequal(base.bus_ids,stress.bus_ids));
limits=[min([base.voltage(:);stress.voltage(:)]),max([base.voltage(:);stress.voltage(:)])];
labels={'基线：负荷倍率1.00，α=1.026','负荷增加10%：α=1.1286'};
data={base,stress};

% Keep an explicit discrete heatmap; each rectangle is one bus-hour sample.
fig=figure('Color','w','Position',[60 60 1150 850]);
tiledlayout(2,2,'TileSpacing','loose','Padding','loose');
for k=1:2
    ax=nexttile;imagesc(1:24,data{k}.bus_ids,data{k}.voltage);axis xy;
    clim(limits);cb=colorbar;cb.Label.String='电压 / p.u.';cb.Label.FontName='SimSun';
    title(labels{k},'FontName','SimSun','FontSize',16);
    xlabel('时段 / h','FontName','SimSun');ylabel('节点编号','FontName','SimSun');
    set(ax,'FontSize',14,'XTick',[1 6 12 18 24],'YTick',[1 10 20 30 40 49]);
end
ax=nexttile;voltage_lines(base,stress);set(ax,'FontSize',14);
ax=nexttile;loss_lines(base,stress);set(ax,'FontSize',14);

end

function voltage_lines(base,stress)
plot(1:24,base.voltage(base.min_bus,:),'-o','LineWidth',1.8,'MarkerSize',5);hold on;
plot(1:24,stress.voltage(stress.min_bus,:),'-s','LineWidth',1.8,'MarkerSize',5);
yline(0.95,'k--','HandleVisibility','off');grid on;xlim([1 24]);
xlabel('时段 / h','FontName','SimSun');ylabel('电压 / p.u.','FontName','SimSun');
title('节点13的24小时电压','FontName','SimSun','FontSize',16);
legend({'基线','负荷增加10%'},'FontName','SimSun','FontSize',13,'Location','best');
end

function loss_lines(base,stress)
plot(1:24,base.loss_mw,'-o','LineWidth',1.8,'MarkerSize',5);hold on;
plot(1:24,stress.loss_mw,'-s','LineWidth',1.8,'MarkerSize',5);grid on;xlim([1 24]);
xlabel('时段 / h','FontName','SimSun');ylabel('有功损耗 / MW','FontName','SimSun');
title('逐时网损：每时段1小时','FontName','SimSun','FontSize',16);
legend({'基线','负荷增加10%'},'FontName','SimSun','FontSize',13,'Location','best');
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