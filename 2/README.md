# 2 — 圆柱绕流多保真（LF / HF）数据集生成（redbKIT）

本目录接在 `../1`（$Re=100$ 单工况基础算例）之后，复现 `../yuanzhuraoliu_LSTM.md`（Conti, Guo, Manzoni, Hesthaven, CMAME 2023）第 5.1 节的多保真数据集：

| 保真度 | 网格 | 时间步长 | Re 取值 | 每个 Re 的输出 |
|---|---|---|---|---|
| **LF** | Gmsh 粗网格 `mesh/LF_coarse.msh` | $\Delta t_{LF}=0.02$ | 19 个，在 $[70,160]$ 内均匀分布：$70,75,\dots,160$ | 26 个时刻 × $(C_D, C_L, F_x, F_y)$ |
| **HF** | Gmsh 细网格 `mesh/HF_fine.msh` | $\Delta t_{HF}=0.01$ | 10 个，在 $[70,160]$ 内均匀分布：$70,80,\dots,160$ | 26 个时刻 × $(C_D, C_L, F_x, F_y)$ |

每个工况都从 $t=0$ 计算到 $t=15$，然后取 $t\in[14.5,15.0]$ 上均匀分布的 26 个时刻 $t_j=14.5+0.02\,j$（$j=0,\dots,25$）作为数据集。两种时间步长下，这 26 个时刻都正好落在时间步上，所以直接取计算值，没有插值。

几何、边界条件、入流曲线、物性、$Re=0.1/\nu$ 的定义以及 P2–P1 有限元 + BDF2 离散，都与 `../1` 和论文完全一致。设置细节见 `../1/README.md` 第 1 节。

## 1. 网格（Gmsh）

`mesh/cylinder.geo` 是参数化的几何文件，与 redbKIT `obstacleDFG.geo` 的几何和物理标记相同：1 入口、2 出口、3 上下壁面、4 圆柱。`s1` 是四个角点处的网格尺寸，`s3` 是圆柱表面的网格尺寸，都通过 `-setnumber` 传入。网格采用 Delaunay 算法（`Mesh.Algorithm = 5`），输出为 msh2.2 格式。

| | `s1` | `s3` | 节点 | 三角形 | 边界线单元 | Gmsh 单元总数 | P2–P1 自由度 | 论文 Fig. 6 |
|---|---|---|---|---|---|---|---|---|
| LF | 0.02 | 0.005 | **3903** | 7480 | 326 | **7806** | 33 253 | 3899 节点 / 7798 单元 |
| HF | 0.015 | 0.002 | **8236** | 15962 | 510 | **16472** | 71 174 | 8239 节点 / 16478 单元 |

两套网格的节点数和单元数与论文 Fig. 6 相差都不到 0.1%。论文的两套网格正是 redbKIT 自带的 `obstacleDFG_Coarse.msh` 和 `obstacleDFG.msh`。见 `results/fig_mesh_LF_HF.png` 和 `results/mesh_info.txt`。

## 2. 配置文件（YAML）

`config/LF.yaml` 和 `config/HF.yaml` 包含一个工况集的全部参数：网格文件与尺寸、物性、有限元、时间离散、Newton 参数、Re 网格、采样窗口、输出目录和检查点间隔。两者只有以下几项不同：

```yaml
# LF.yaml                         # HF.yaml
mesh:  file: mesh/LF_coarse.msh   # mesh/HF_fine.msh
       s1: 0.02, s3: 0.005        # s1: 0.015, s3: 0.002
time:  dt: 0.02                   # dt: 0.01
reynolds: {min: 70, max: 160, n: 19}   # n: 10
checkpoint: every_steps: 25       # 50   （都等于每 0.5 个时间单位存一次）
```
YAML 由 `src/read_yaml.m` 读取（纯 MATLAB/Octave 实现，不需要工具箱）。

## 3. 批量求解与断点续算

```
run_batch('config/LF.yaml')                      % 19 个 LF 工况
run_batch('config/HF.yaml')                      % 10 个 HF 工况
run_batch({'config/HF.yaml','config/LF.yaml'})   % 全部 29 个工况（本次实际使用）
run_batch('config/HF.yaml', [100 110])           % 只算指定的 Re
```

`run_batch` 依次对每个 Re 调用 `src/solve_case.m`。每个工况的输出放在 `results/<LF|HF>/Re<Re>/` 下：

| 文件 | 内容 |
|---|---|
| `forces.txt` | 完整时间历程 $t, C_D, C_L$，计算过程中持续写入 |
| `checkpoint.mat` | 求解器状态：当前速度/压力，BDF2 所需的历史解，已算出的力历程，Newton 迭代次数。每 `every_steps` 步原子写入一次（先写 `*_tmp.mat` 再 `movefile`）。工况完成后删除 |
| `result.mat` | 结构体 `R`，包括 `t, CD, CL, Fx, Fy`（全历程）、`t_s, CD_s, CL_s, Fx_s, Fy_s`（26 个采样点）、`newton_its`、网格信息和耗时。**它的存在就表示该工况已完成** |
| `sessions.txt` | 每次启动/续算的记录（续算起点、耗时） |

**断点续算的机制：**
1. 已有 `result.mat` 的工况直接跳过，不会重算。
2. 中断的工况（崩溃、卡死、Ctrl‑C、关机等）从最后一个 `checkpoint.mat` 继续。`src/NSt_Solver_ckpt.m` 是 redbKIT `NSt_Solver.m` 的修改版，会恢复 BDF 的历史状态。实测中，在中途 `kill -9` 后再续算，得到的力系数与不中断的计算**逐位相同**（最大差 0）。
3. 恢复的方法就是重新执行同一条命令。

**并行 worker：** 多个 MATLAB/Octave 进程可以同时运行 `run_batch`。正在计算的工况会写一个 `lock.txt`（记录进程号），其他进程自动跳过它，去计算下一个。进程已退出、或 20 分钟没有进展的锁视为失效。

**看门狗（Linux）：** `./run_all.sh [N_WORKERS]` 依次完成以下步骤：
1. 生成网格；
2. 启动 N 个单线程 Octave worker；
3. 每分钟检查一次：worker 崩溃就重启；日志 20 分钟没有更新（判定为卡死）就 `kill -9` 并从检查点重启；
4. 29 个工况全部完成后，自动执行 `build_datasets` 和 `plot_results`。

本次运行：2 个 worker，2 核 CPU，Octave 8.4。

| 项目 | 数值 |
|---|---|
| 运行时段 | 2026‑11‑02 10:34 → 23:50 |
| 墙钟时间 | **13.3 h** |
| 累计计算 | 26.2 CPU·h |
| HF 平均耗时 | 114 min / 工况 |
| LF 平均耗时 | 22.6 min / 工况 |
| Newton 迭代 | 每步 3–4 次，全部收敛（tol $10^{-8}$） |

日志在 `results/run_logs/`。

## 4. 数据集

运行 `build_datasets` 生成：

| 文件 | 内容 |
|---|---|
| `results/dataset_LF.mat` | `Re` (19×1)，`t` (1×26)，`CD`、`CL`、`Fx`、`Fy` (19×26)，以及 `dt, mesh, nodes, triangles, dofs, description` |
| `results/dataset_HF.mat` | 同上，Re 维度为 10 |
| `results/dataset_LF.csv` / `dataset_HF.csv` | 长表格式，每行一个 $(Re,t)$，列为 `Re,t,CD,CL,Fx,Fy`（494 / 260 行数据） |

定义：$C_D=\dfrac{2F_x}{\rho U_{mean}^2D}$，$C_L=\dfrac{2F_y}{\rho U_{mean}^2D}$，其中 $\rho=1$、$U_{mean}=1$、$D=0.1$，因此 $F_{x,y}=0.05\,C_{D,L}$。$F_x$、$F_y$ 是作用在圆柱上的总力（含压力和黏性力）。$C_L$ 保留正负号；论文中的 $|C_L|$ 就是 `abs(CL)`。读取方法：MATLAB/Octave 中用 `D = load('results/dataset_HF.mat')`，Python 中用 `pandas.read_csv(...)` 或 `scipy.io.loadmat(...)`。

## 5. 结果

### 5.1 LF 与 HF 在同一 Re 下的对比（$Re=100$）
`results/fig_CD_CL_LF_vs_HF_Re100.png` 对比 $t\in[14.5,15]$ 上的结果（曲线是全部时间步，符号是 26 个数据集采样点）。`results/fig_history_LF_vs_HF_Re100.png` 是 $t\in[0,15]$ 的全历程。

| $Re=100$，$t\in[10,15]$ | LF | HF | LF 相对 HF | `../1`（redbKIT 原网格） |
|---|---|---|---|---|
| $\overline{C_D}$ | 3.1523 | 3.1889 | −1.1 % | 3.1897 |
| $C_{D,\max}$ | 3.1757 | 3.2196 | −1.4 % | 3.2195 |
| $C_{L,\max}$ | 0.8239 | 0.9600 | −14.2 % | 0.9600 |
| $St=fD/U_{mean}$ | 0.2885 | 0.2979 | −3.2 % | 0.2979 |

HF 结果与 `../1` 基本一致（$C_{D,\max}$ 相差 $10^{-4}$），说明 Gmsh 细网格可以代替 redbKIT 原网格。LF 由于网格较粗、时间步较大，有系统偏差：数值耗散偏大，因此升力幅值偏小、St 偏小、阻力略低。到 $t=14.5$ 时，相位误差已累积到接近半个周期。这正是多保真模型要学习的 LF→HF 修正。

### 5.2 全部 Re
- `results/fig_LF_HF_maps.png`：仿照论文 Fig. 7，给出 LF（19 Re）和 HF（10 Re）的 $C_L(t,Re)$ 与 $C_D(t,Re)$ 分布图，色标相同。
- `results/fig_LF_vs_HF_allRe.png`：在 10 个共同 Re 上，LF 与 HF 的 $C_D$、$C_L$ 曲线。
- `results/fig_summary_vs_Re.png` 和 `results/summary.csv`：$\overline{C_D}$、$C_{L,\max}$、St 随 Re 的变化。

| Re | 70 | 80 | 90 | 100 | 110 | 120 | 130 | 140 | 150 | 160 |
|---|---|---|---|---|---|---|---|---|---|---|
| HF $\overline{C_D}$ | 3.391 | 3.300 | 3.235 | 3.189 | 3.158 | 3.138 | 3.126 | 3.121 | 3.118 | 3.118 |
| LF $\overline{C_D}$ | 3.373 | 3.278 | 3.207 | 3.152 | 3.112 | 3.081 | 3.060 | 3.046 | 3.036 | 3.029 |
| HF $C_{L,\max}$ | 0.553 | 0.692 | 0.825 | 0.960 | 1.099 | 1.240 | 1.376 | 1.505 | 1.623 | 1.731 |
| LF $C_{L,\max}$ | 0.489 | 0.612 | 0.720 | 0.824 | 0.930 | 1.036 | 1.137 | 1.249 | 1.352 | 1.447 |
| HF St | 0.2848 | 0.2905 | 0.2948 | 0.2979 | 0.3003 | 0.3024 | 0.3043 | 0.3061 | 0.3079 | 0.3097 |
| LF St | 0.2759 | 0.2813 | 0.2854 | 0.2885 | 0.2908 | 0.2925 | 0.2939 | 0.2952 | 0.2963 | 0.2975 |

所有工况在 $t\ge 10$ 后都达到稳定的周期性涡脱落。随 Re 增大，$\overline{C_D}$ 单调减小，$C_{L,\max}$ 和 St 单调增大；LF 在整个参数区间内的偏差平滑且有规律。

## 6. 文件说明

| 文件 | 作用 |
|---|---|
| `config/LF.yaml`, `config/HF.yaml` | 两套计算配置 |
| `mesh/cylinder.geo` | 参数化 Gmsh 几何 |
| `mesh/LF_coarse.msh`, `mesh/HF_fine.msh` | 生成的粗、细网格 |
| `generate_meshes.m` | 读取 YAML，调用 `gmsh` 生成两套网格，输出统计信息和网格图 |
| `run_batch.m` | 批量求解，支持断点续算和多进程锁 |
| `run_all.sh` | Linux 一键运行 + 看门狗 |
| `build_datasets.m` | 汇总成 `dataset_LF/HF.mat/.csv` |
| `plot_results.m` | 所有对比图和 `summary.csv` |
| `src/solve_case.m` | 单个 (保真度, Re) 工况：求解、续算、采样、保存 |
| `src/NSt_Solver_ckpt.m` | redbKIT `NSt_Solver` + 检查点/续算 + 简洁日志 + 返回力历程 |
| `src/datafile_cylinder.m` | redbKIT 数据文件，参数来自 YAML |
| `src/read_yaml.m` | 轻量 YAML 读取器 |
| `src/setup_paths.m` | 添加 redbKIT 路径（环境变量 `REDBKIT_PATH`，默认 `../../redbKIT`） |
| `src/load_msh_2D.m`, `src/octave_compat/` | 与 `../1` 相同：msh 读取器和 Octave 兼容函数 |
| `results/` | 数据集、图、`summary.csv`、`mesh_info.txt`、各工况结果、`run_logs/` |

## 7. 复现

**MATLAB：**
```matlab
% 先下载 redbKIT 并在其根目录执行一次 >> make；Gmsh >= 4 需在系统 PATH 中
setenv('REDBKIT_PATH', '/path/to/redbKIT');
cd 2
generate_meshes                                 % 已提供 .msh，可以跳过
run_batch({'config/HF.yaml','config/LF.yaml'})  % 中断后重新执行即可续算
build_datasets
plot_results
```
如果要在 MATLAB 中并行计算，可以同时打开多个 MATLAB 窗口执行同一条 `run_batch` 命令。

**GNU Octave（本次实际运行方式）：** 按 `../1/setup_and_run_octave.sh` 编译好 redbKIT 后，执行
```bash
REDBKIT_PATH=/path/to/redbKIT ./run_all.sh 2
```

> 注意：已经完成的工况会被跳过。如果修改了配置后需要重算，请先删除 `results/<LF|HF>/` 中对应的 `Re*/` 目录。
