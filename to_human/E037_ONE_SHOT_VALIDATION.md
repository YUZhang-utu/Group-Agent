# E037 一次执行后续验证

台式机可用时，直接运行本套件，无需先单独运行 E036 或重跑百万 E035。

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e037_workstation_suite.sh
```

套件默认使用原 E034 目录和原大库路径，输出：

```text
/mnt/local/hand/yuzhang/aidd/e037-workstation-suite/run-v1
```

## 会运行什么

| 配置 | coarse / refine chunk | pair-seed 生成 |
|---|---|---|
| old-reference | 2000 / 250 | 原实现：全部生成再截取 |
| fine-reference | 500 / 64 | 原实现 |
| fine-bounded | 500 / 64 | 达到512个唯一seed时停止 |
| old-bounded | 2000 / 250 | 达到512个唯一seed时停止 |

每个配置默认16workers，每次执行两条已校准 WEE1 查询。四组运行两轮，第二轮反序，
共16次查询执行。另有一组 fine-bounded profiling（两条查询），只分析 coarse/refine
第一块，计时不纳入速度比较。第一块只是诊断样本，不能代表所有候选的热点分布。
这不是减少 seed 数：所有配置都保留相同的前512个唯一seed，顺序也不变。

每组都从已建好的全库 FAISS 索引检索候选，再做 Gaussian 和优化版 E031。
预算10000、Top-N5000、nprobe128、评分方程和 E031 注释用途保持不变。
检索候选必须精确匹配 E034；Gaussian 的每个输出数组、变换矩阵和 E031 分数/
分配/排名都要通过等价性检查。失败会记录日志、继续其他组，最终套件返回非零状态。

## 输出与恢复

最终发回根目录 `report.md` 和 `report.json`。汇总按查询列出各配置的执行时间中位数、
相对 old-reference 的速度比和计时有效性。每轮顺序已反转，但样本量仍小、系统负载
无法控制，属于探索性工程测量，不能视为服务 SLA。

`progress.json` 记录已完成组数；根目录每组有对应 `.log` 文件，可在另一个终端
用 `tail -f` 查看当前文件。子目录中有单次报告和原始输出。
profiling 文件在 `profile-fine-bounded/{8bju,1x8b}/gaussian/{coarse,refine}/` 下：
`worker-first-chunk-profile.txt` 和 `worker-first-chunk.pstats`。

中断后用同一命令加 `--resume`：

```bash
bash scripts/run_e037_workstation_suite.sh --resume
```

完成的查询保留原计时，未完成的 chunk 可恢复；发生部分复用的查询不参与速度比。
要获得完整的新计算计时，可指定新的 `E037_OUTPUT` 再运行。
不要更新代码后向旧协议目录混入新结果；代码/参数/输入变化会被拒绝。
不需要删除或覆盖 E034、E035 或已有 E036 输出。

环境变量：`E037_WORKERS`（16）、`E037_REPEATS`（2，至少2）、`E037_OUTPUT`、
`E037_E034_SOURCE`、`AIDD_BATCH`。套件期间不要同时运行其他重计算实验。
启动完整性检查与索引加载单独记录，原始秒数在子报告中；不会运行全库精确参考扫描。

## 本机已完成

- 205测试通过，2项缺失依赖跳过；包括多进程、恢复、损坏拒绝、套件失败后继续、
  逆序调度、seed精确前缀、完整Gaussian数组一致和独立profiling。
- 固定随机种子的12组合成seed生成基准、每组3次：全部前缀严格一致。
  生成量超过512的组，生成器加速1.27–9.69倍；未超过时约0.99–1.01倍。
  原始数据：`to_human/E037_LOCAL_SEED_BENCHMARK.json`。
- 以上仅说明生成器等价且在部分合成场景节省工作，不能推断整条搜索提速幅度。
  真实库的收益、负载平衡和热点需要这次台式机批量结果确认。

原 E035 导出的39/41个分歧姿态仍待人工审查，本套件不替代姿态质量或活性富集验证。
