%% 任务二：储能套利、互斥约束及功率容量比较
% 储能231 雷宇杰 206231313。2 MW/2 MWh，无光伏、无网络约束。 采用LP/MILP及独立矩阵模型，比较2/2、3/2、2/3三个场景。

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

Task2_Result=task2_experiment(task_dir);
disp(Task2_Result.summary);
%% 手算：购入1 MWh用于充电，可出售多少电量

E_stored_hand=1*0.90;
E_sale_hand=E_stored_hand*0.95;
disp(table(E_stored_hand,E_sale_hand));
assert(abs(E_sale_hand-0.855)<1e-12);
fprintf('TASK2_COMPLETE: 12 LP/MILP solves.\n');
function result=task2_experiment(out_dir)
%TASK2_EXPERIMENT Verify three personal arbitrage scenarios via two formulations.
assert(startsWith(which('linprog'),matlabroot,'IgnoreCase',true));
assert(startsWith(which('intlinprog'),matlabroot,'IgnoreCase',true));
profiles=load('profiles_24h.mat');
assert(profiles.T==24 && all(profiles.price_MWh>0));
names={'Base_P2_E2','Power_P3_E2','Energy_P2_E3'};
powers=[2 3 2];energies=[2 2 3];
method_names={'YALMIP_LP','YALMIP_MILP','Matrix_LP','Matrix_MILP'};
rows=cell(12,14);times=cell(12,1);n=0;
scenarios=struct([]);
for k=1:3
    p=struct('price',profiles.price_MWh(:),'power',powers(k),'energy',energies(k), ...
        'eta_ch',0.90,'eta_dis',0.95,'soc0',0.50,'smin',0.10,'smax',1.00,'dt',1);
    scenario=struct('name',names{k},'parameters',p,'solutions',struct([]));
    for method=1:4
        if method<=2
            r=solve_yalmip(p,method==2);
        else
            r=solve_matrix(p,method==4);
        end
        r.method=method_names{method};
        r.check=verify_task2(r,p);
        if method==1
            scenario.solutions=r;
        else
            scenario.solutions(method)=r;
        end
        n=n+1;
        rows(n,:)={string(names{k}),string(r.method),p.power,p.energy,string(r.solver), ...
            r.revenue_yuan,sum(r.pch_mw)*p.dt,sum(r.pdis_mw)*p.dt, ...
            r.check.energy_balance_error,r.check.soc_terminal_error, ...
            r.check.simultaneous_mw,min(r.soc),max(r.soc),r.seconds};
        times{n}=table(repmat(string(names{k}),24,1),repmat(string(r.method),24,1), ...
            (1:24)',p.price,r.pch_mw,r.pdis_mw,r.soc, ...
            'VariableNames',{'Scenario','Method','Hour','Price_YuanPerMWh','Charge_MW','Discharge_MW','SoC'});
        fprintf('TASK2 %-14s %-12s revenue=%.9f energy_res=%.3e terminal=%.3e overlap=%.3e time=%.3fs\n', ...
            names{k},r.method,r.revenue_yuan,r.check.energy_balance_error, ...
            r.check.soc_terminal_error,r.check.simultaneous_mw,r.seconds);
    end
    objectives=[scenario.solutions.revenue_yuan];
    assert(max(objectives)-min(objectives)<1e-5, ...
        'Task2:IndependentAgreement','LP/MILP and independent matrix objectives disagree.');
    scenario.max_revenue_difference=max(objectives)-min(objectives);
    if k==1
        scenarios=scenario;
    else
        scenarios(k)=scenario; %#ok<AGROW>
    end
end
base=scenarios(1).solutions(2);power3=scenarios(2).solutions(2);energy3=scenarios(3).solutions(2);
assert(power3.revenue_yuan>=base.revenue_yuan-1e-5 ...
    && energy3.revenue_yuan>=base.revenue_yuan-1e-5, ...
    'Task2:Monotonic','Relaxed power/energy limits reduced optimal revenue.');
summary=cell2table(rows,'VariableNames',{'Scenario','Method','Power_MW','Energy_MWh', ...
    'Solver','Revenue_Yuan','Charge_MWh','Discharge_MWh','EnergyBalanceError_MWh', ...
    'TerminalSoCError','SimultaneousPower_MW','MinSoC','MaxSoC','SolveSeconds'});
checks=[scenarios.solutions];checks=[checks.check];
metrics=struct('student_id',206231313,'base_revenue_yuan',base.revenue_yuan, ...
    'lp_revenue_yuan',scenarios(1).solutions(1).revenue_yuan, ...
    'milp_revenue_yuan',base.revenue_yuan,'power3_revenue_yuan',power3.revenue_yuan, ...
    'energy3_revenue_yuan',energy3.revenue_yuan, ...
    'energy_balance_error',max([checks.energy_balance_error]), ...
    'soc_terminal_error',max([checks.soc_terminal_error]), ...
    'simultaneous_mw',max([checks.simultaneous_mw]), ...
    'matrix_lp_revenue_yuan',scenarios(1).solutions(3).revenue_yuan, ...
    'matrix_milp_revenue_yuan',scenarios(1).solutions(4).revenue_yuan, ...
    'max_four_solver_revenue_difference',max([scenarios.max_revenue_difference]), ...
    'base_charge_mwh',sum(base.pch_mw),'base_discharge_mwh',sum(base.pdis_mw), ...
    'base_min_soc',min(base.soc),'base_max_soc',max(base.soc), ...
    'yalmip_version',yalmip('version'),'matlab_version',version, ...
    'native_linprog',which('linprog'),'native_intlinprog',which('intlinprog'), ...
    'generated_at',char(datetime('now')));
result=struct('student_id',206231313,'name','雷宇杰','class','储能231', ...
    'metrics',metrics,'scenarios',scenarios,'summary',summary);
save_week1_part(out_dir,'Task2_Result',result);
plot_baseline(scenarios(1),out_dir);
fprintf('TASK2 revenues: base=%.9f power3=%.9f energy3=%.9f yuan\n', ...
    base.revenue_yuan,power3.revenue_yuan,energy3.revenue_yuan);
end

function r=solve_yalmip(p,is_milp)
yalmip('clear');T=numel(p.price);
pc=sdpvar(T,1);pd=sdpvar(T,1);s=sdpvar(T,1);
c=[0<=pc<=p.power,0<=pd<=p.power,p.smin<=s<=p.smax,s(T)==p.soc0];
previous=[p.soc0;s(1:T-1)];
c=[c,p.energy*(s-previous)==(p.eta_ch*pc-pd/p.eta_dis)*p.dt];
if is_milp
    u=binvar(T,1);c=[c,pc<=p.power*u,pd<=p.power*(1-u)];
    solver='intlinprog';
else
    solver='linprog';
end
objective=sum(p.price.*(pc-pd))*p.dt;
ops=sdpsettings('solver',solver,'verbose',0,'savesolveroutput',1);
if is_milp
    ops.intlinprog.RelativeGapTolerance=1e-10;
    ops.intlinprog.AbsoluteGapTolerance=1e-8;
else
    ops.linprog.ConstraintTolerance=1e-9;
    ops.linprog.OptimalityTolerance=1e-9;
end
t=tic;sol=optimize(c,objective,ops);seconds=toc(t);
assert(sol.problem==0,'Task2:YALMIP','%s failed: %s',solver,sol.info);
if is_milp,binary=value(u);else,binary=[];end
r=make_result(value(pc),value(pd),value(s),binary,-value(objective), ...
    is_milp,solver,sol.problem==0,sol.info,seconds);
end

function r=solve_matrix(p,is_milp)
% Independent formulation eliminates SoC variables using cumulative energy.
% x=[Pch(1:T);Pdis(1:T)] (and binary u for MILP), all grid-side MW.
T=numel(p.price);C=tril(ones(T));
E=p.dt*[p.eta_ch*C -C/p.eta_dis];
A=[E;-E];b=[repmat(p.energy*(p.smax-p.soc0),T,1);repmat(p.energy*(p.soc0-p.smin),T,1)];
Aeq=E(end,:);beq=0;
f=p.dt*[p.price;-p.price];lb=zeros(2*T,1);ub=repmat(p.power,2*T,1);
if is_milp
    A=[A zeros(2*T,T);eye(T) zeros(T) -p.power*eye(T); ...
        zeros(T) eye(T) p.power*eye(T)];
    b=[b;zeros(T,1);repmat(p.power,T,1)];
    Aeq=[Aeq zeros(1,T)];f=[f;zeros(T,1)];
    lb=[lb;zeros(T,1)];ub=[ub;ones(T,1)];
    options=optimoptions('intlinprog','Display','off', ...
        'RelativeGapTolerance',1e-10,'AbsoluteGapTolerance',1e-8);
    t=tic;[x,fvalue,flag,output]=intlinprog(f,2*T+(1:T),A,b,Aeq,beq,lb,ub,options);seconds=toc(t);
    solver='intlinprog';
else
    options=optimoptions('linprog','Display','none', ...
        'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
    t=tic;[x,fvalue,flag,output]=linprog(f,A,b,Aeq,beq,lb,ub,options);seconds=toc(t);
    solver='linprog';
end
assert(flag==1,'Task2:Native','%s independent matrix solve failed: %s',solver,output.message);
pc=x(1:T);pd=x(T+(1:T));s=p.soc0+E*x(1:2*T)/p.energy;
if is_milp,binary=x(2*T+(1:T));else,binary=[];end
r=make_result(pc,pd,s,binary,-fvalue,is_milp,solver,flag==1,output.message,seconds);
end

function r=make_result(pc,pd,s,binary,revenue,is_milp,solver,success,info,seconds)
r=struct('pch_mw',pc(:),'pdis_mw',pd(:),'soc',s(:),'binary_state',binary(:), ...
    'revenue_yuan',revenue,'is_milp',is_milp,'solver',solver, ...
    'solver_success',success,'solver_info',info,'seconds',seconds);
end

function plot_baseline(scenario,out_dir)
p=scenario.parameters;r=scenario.solutions(2);
fig=figure('Color','w','Position',[60 60 1080 800]);
try,theme(fig,'light');catch,end
tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
nexttile;stairs(1:24,p.price/1000,'k-','LineWidth',1.8);grid on;xlim([0.5 24.5]);ylim([0.25 1.25]);
ylabel('电价 / 元·kWh^{-1}','FontName','SimSun');title('江苏分时电价（课程输入）','FontName','SimSun');
nexttile;bar(1:24,[r.pch_mw -r.pdis_mw],'grouped');grid on;xlim([0.5 24.5]);
ylabel('端口功率 / MW','FontName','SimSun');
legend({'充电（正）','放电（负）'},'FontName','SimSun','Location','best');
title(sprintf('MILP基线：2 MW / 2 MWh；日净电费收益 %.2f 元',r.revenue_yuan),'FontName','SimSun');
nexttile;plot(0:24,100*[p.soc0;r.soc],'-o','LineWidth',1.8,'MarkerSize',4);hold on;
yline(10,'k--');yline(100,'k--');yline(50,':','首末50%','FontName','SimSun');
ylim([0 110]);xlim([0 24]);grid on;
xlabel('时段末 / h','FontName','SimSun');ylabel('SoC / %','FontName','SimSun');
title('充电效率0.90；放电效率0.95；初末SoC=0.50','FontName','SimSun');
sgtitle('任务二  雷宇杰 206231313 | 无光伏、无网络约束','FontName','SimSun');
set(findall(fig,'Type','axes'),'FontName','Times New Roman','FontSize',16);
set(findall(fig,'Type','legend'),'FontSize',14);
ess_live_export(fig,out_dir,'Fig_Week1_Task2_206231313');
end

function check=verify_task2(r,p)
% Recompute physics and economics from numerical schedules, not solver status.
pc=r.pch_mw(:);pd=r.pdis_mw(:);s=r.soc(:);T=numel(p.price);
assert(numel(pc)==T && numel(pd)==T && numel(s)==T ...
    && all(isfinite([pc;pd;s])),'Task2:Shape','Invalid schedule.');
assert(r.solver_success,'Task2:Solver','A solver did not certify its LP/MILP solve.');
previous=[p.soc0;s(1:end-1)];
residual=p.energy*(s-previous) ...
    -(p.eta_ch*pc-pd/p.eta_dis)*p.dt;
check.energy_balance_error=max(abs(residual));
check.soc_terminal_error=abs(s(end)-p.soc0);
check.simultaneous_mw=max([0;min(pc,pd)]);
check.power_violation=max([0;-pc;-pd;pc-p.power;pd-p.power]);
check.soc_violation=max([0;p.smin-s;s-p.smax]);
check.revenue_yuan=sum(p.price(:).*(pd-pc))*p.dt;
assert(check.energy_balance_error<1e-7,'Task2:Energy','SoC energy balance failed.');
assert(check.soc_terminal_error<1e-8,'Task2:Terminal','Terminal SoC failed.');
assert(check.power_violation<1e-7 && check.soc_violation<1e-7, ...
    'Task2:Bounds','Power/energy bounds failed.');
assert(abs(check.revenue_yuan-r.revenue_yuan)<1e-6, ...
    'Task2:Revenue','The saved revenue has the wrong sign or units.');
assert(check.simultaneous_mw<1e-6,'Task2:Overlap', ...
    'Positive-price LP or explicit MILP has simultaneous charging/discharging.');
assert(abs(sum(pd)-p.eta_ch*p.eta_dis*sum(pc))<1e-7, ...
    'Task2:Cycle','Cyclic grid-side energy accounting failed.');
if r.is_milp
    assert(numel(r.binary_state)==T ...
        && max(abs(r.binary_state-round(r.binary_state)))<1e-6, ...
        'Task2:Integer','Binary state is not integral.');
    assert(max(pc-p.power*r.binary_state)<1e-7 ...
        && max(pd-p.power*(1-r.binary_state))<1e-7, ...
        'Task2:BinaryPower','Charging/discharging status disagrees with power.');
end
end

function ess_live_export(fig,folder,name)
% 按课程要求导出白底 PNG、矢量 PDF 和 FIG。
if ~isfolder(folder), mkdir(folder); end
exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',600,'BackgroundColor','white');
exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector','BackgroundColor','white');
savefig(fig,fullfile(folder,[name '.fig']));
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