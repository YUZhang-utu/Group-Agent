# E034：扩大库 WEE1 检索 → Gaussian 精修 → E031

代码用于下一阶段实测，尚不代表 25.8M 库上的 E034 已通过。
E033 已收到的结果为 25,813,808 构象、8,318,351 个来源分组分子，
8 查询校准面板通过；10,000 候选 / nprobe=128 是暂定起点。

本地验证：166 项通过，2 项依赖相关跳过；其中新增19项覆盖 E034，
包含真实 mmap 文件上的 Gaussian→E031 数值衔接。Bash 语法检查通过。
当前 Windows 解释器缺少 RDKit/FAISS，真实化学查询准备及工作站全流程
仍需在原 Linux 环境运行；不将模拟索引测试当成真实 FAISS 验收。

## 台式机执行

在之前运行 E032/E033 的 Python 环境中：

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e034_expanded_wee1.sh
```

默认输入：

- E032：`/mnt/local/hand/yuzhang/aidd/library-precompute-20260911`
- E033：`/mnt/local/hand/yuzhang/aidd/e033-library-acceptance/20260918-091604`
- QT9：仓库 `data/e019_query_8bju/{8BJU.cif,QT9.cif,query_manifest.json}`
- 824：仓库 `data/e026_query_1x8b/{1X8B.cif,824.cif,query_manifest.json}`
- 默认输出：`/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/20260918-run1`

查询是之前 E026 已使用的两个本地查询目录；数据不在 Git 中。
如果目录不同，只指定真实路径，不能把另一条查询的文件改名替代：

```bash
E034_QT9_DIR=/actual/qt9-query-dir \
E034_1X8B_DIR=/actual/1x8b-query-dir \
bash scripts/run_e034_expanded_wee1.sh
```

脚本不联网下载结构，也不重建库。缺少原始查询文件、库未完成、E033
报告损坏/校准未过、查询身份不符或库来源变更时明确报错。
依赖与原 E032/E033/E026 环境相同：Python 3.11+、NumPy、RDKit、
Biopython、Gemmi、FAISS。默认 Gaussian 16 个 worker、FAISS 20 线程，
BLAS 单线程；可用 E034_WORKERS/E034_THREADS 覆盖。
大量分片 mmap 需要足够的文件描述符。脚本会在既有硬上限内提高当前
进程的软上限（子进程继承），不修改系统全局配置；硬上限不足则提前报错。

## 执行内容与边界

1. 读取 E033 完成标记并校验 report.json；核对协议、库清单、变换、索引
   哈希、当前登记/构象/化学伴随数据 ID。当前全量 payload 字节验收继承
   已完成 E033，本次重查元数据/身份及索引字节，不重复读取全部原始 MOL2。
2. 保存锁定查询输入副本；从真实晶体坐标生成 USRCAT 和 Gaussian 查询。
   USRCAT 使用已有标准化 parent 路径，Gaussian 保留 CCD 配体化学形式，
   两者来自同一晶体配体实例。外部查询无已验证库 molecule-ID 映射，
   不声称已经按化学结构去除全库自匹配。
3. 一次全库扫描生成两条查询的精确 USRCAT Top-1000，分别测 10,000
   候选下 nprobe=128/256。两条查询严格 Top-1000 召回均 >=95% 时选128；
   否则尝试256；均失败则只生成检索报告，退出码2，不继续精修。
4. 对选定档位的完整构象候选做 Gaussian 粗评分和每目标 Top-5000
   并集精修，保持原有512 pair-seed上限。不会提前每分子只留一个构象。
5. 在原有三套刚体姿态上运行 E031 侧评分；保存相互作用覆盖、逐 anchor
   匹配、相关性和 Top-K 重合度，不改原有排名、不做 docking。
6. 分别报告检索、候选读取、精确参考、Gaussian、E031 耗时和保留数量。
   E031 <=10% 门槛用同次同查询的 Gaussian 精修阶段时间。未通过也保存
   结果，执行 complete 不等于延迟门槛或科学结论通过。

搜索参数的靶点校准不等于 Gaussian 或活性召回；预算缩减不等于化学淘汰。
未引入药效团额外候选通道，以便单独度量校准后的 USRCAT→精修路径。
旧 `run_e031_key_interaction_matching.sh` 不再用固定77.99秒输出门槛判断。

## 中断续跑

```bash
bash scripts/run_e034_expanded_wee1.sh --resume
```

同一输出目录必须匹配输入、源码、依赖版本、机器和 worker/thread 参数。
已完成阶段检查哈希后复用，保留原始计时；中断的 Gaussian 阶段复用其
已验证 chunk。若精修阶段只有部分重新计算，本轮延迟门槛标记 unavailable，
不能把剩余计算时间当成完整基线。获取新的完整计时使用新输出：

```bash
E034_OUTPUT=/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/timing-run2 \
bash scripts/run_e034_expanded_wee1.sh
```

原有输入、结果不会自动删除。阶段间掉电可能需要重做尚未写出完成回执的
短阶段。FAILED.json 保留历史失败原因；当前执行状态看 RUN_STATUS.json，
并与最终完成标记及 report.json 状态核对。

## 返回哪些内容

优先返回输出目录中的 `report.md`、`report.json`。
完整候选和精确距离保存在 `retrieval/`；每条查询下有输入副本、
`gaussian/refine/merged-scores.npz`、`interaction-matches.npz` 及 manifest。
根目录 `protocol.json` 保存输入/代码指纹，`*.stage.json` 保存阶段回执。
索引加载、验收与精确参考扫描耗时不计入在线查询速度；小样本计时不是 SLA。
