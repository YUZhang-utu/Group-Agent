# E033：登记后验收、速度／质量／压缩比例评测

## 台式机执行

在原来可以运行 E032 的 Python 环境中：

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e033_library_acceptance.sh
```

脚本默认读取：
`/mnt/local/hand/yuzhang/aidd/library-precompute-20260911`。
报告写入独立的、带时间戳的同级 `e033-library-acceptance-*` 目录。
不重建索引，不修改 MOL2、登记数据库或预计算产物，也不会上传数据。

**登记完成不等于 E032 全部完成。** 若提示缺少 `COMPLETE.json`，先确认原进程
是否仍在生成特征、训练或合并 FAISS；还在运行就让它完成。若原进程已退出，
按下面原命令续跑，完成后再运行上述评测脚本。不要同时启动两个预计算进程。

```bash
set -o pipefail
python scripts/precompute_library_batch.py \
  --source /mnt/local/hand/yuzhang/aidd/mc_data \
  --old-artifacts "$PWD/data/e019_artifacts/catalog.json" \
  --old-chemical /mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json \
  --output /mnt/local/hand/yuzhang/aidd/library-precompute-20260911 \
  --workers 20 --preserve-index-conflicts --run \
  2>&1 | tee -a /mnt/local/hand/yuzhang/aidd/library-precompute-20260911/run.log
```

## 实际评测内容

1. **库验收**：完成标记、源文件登记覆盖、记录数=插入数+精确重复数、登记与
   各阶段构象数量、连续全局 ID、登记／构象数据／化学伴随数据的 ID 一致性、
   清单来源、生成文件与药效团 posting 的完整 SHA256、最终索引数量和维度。
   旧产物保存的原路径通过唯一文件 stem 与来源哈希对应，不要求重写旧清单。
   原始 MOL2 不再全部读一遍；报告明确说明没有重新计算其当前字节哈希。
2. **速度**：索引加载、FAISS 单查询首轮和预热后的 p50/p95、候选读取与精确
   USRCAT 重排耗时、精确参考扫描总耗时、进程峰值 RSS。查询延迟不包含一次性
   全库验收和参考答案生成。首次调用不保证是冷缓存，小样本 p95 不是服务 SLA。
3. **检索质量**：用冻结标准化参数，对完整库分块精确计算 USRCAT 平方 L2，
   比较 FAISS 候选对真实 Top-100、Top-1000 的覆盖。距离相同的边界并列项单独
   报告。不会把某个分片的精确结果当成全库真值。
4. **压缩比例**：逐查询报告保留构象数、独立分子 ID 数、相对全库的压缩比例；
   比较 1,000 / 10,000 / 100,000 个候选预算与 nprobe 64 / 128 / 256。
   同时给出保留 100 / 1,000 / 10,000 个分子的容量情景，每分子一个描述符排名
   最佳的代表构象，并报告对精确邻居分子集合的覆盖。**这些是容量情景，不是
   已经过科学验收的 docking 入选名单，也不是最终应只留一个姿态的结论。**

默认从均匀抽样的构象 ID 中确定性选出 8 个不同分子作为查询（构象较多的分子
被抽中的概率更高，并非分子等概率抽样），并从精确真值与搜索结果中排除查询分子
的全部构象，避免“找回自己”造成虚高召回。这仍然是库内校准集，不是独立的
生物学测试集。工程门槛暂设为所有查询 Top-1000 严格召回 >=95%，先挑通过门槛
的最小候选预算，再比较速度；未通过就明确报告没有合格配置。

## 如何读输出

- `report.md`：优先看这个，列出速度、最差查询召回、剩余分子数和压缩范围。
- `report.json`：完整参数、来源标识、逐查询结果、临时配置建议与能力边界。
- `metrics.csv`：比较每个查询／配置／分子容量情景。
- `acceptance.json`：数据完整性验收。
- `protocol.json`、`queries.npz`、`truth_*.npz`：冻结参数、查询和精确参考答案。
- `candidates_*.npz`：候选全局 ID、分子 ID、精确描述符距离与各容量情景的代表。
  使用未截断的 `global_ids` 可对接既有 Gaussian 候选读取接口；必须搭配同一个
  查询和同一个库，不能把库内抽样查询的候选配给另一条 WEE1 共晶查询。
- `EVALUATION_COMPLETE.json`：评测执行完成；注意这不等于召回门槛通过，查看
  `calibration_status`。`FAILED.json` 表示执行中断，保留已输出文件用于诊断。

先把 `report.md` 和 `report.json` 的内容或文件发回来。不要只发一个平均过滤比例；
我们要同时看最差查询的召回以及剩余的独立分子数量。

## 不能由本次报告直接得出的结论

这一轮评测的是 **3D 描述符粗检索**。USRCAT 近邻一致性不等于 Gaussian 姿态
匹配、关键相互作用保留、活性富集或 docking 成功率。预算导致的候选缩减也不
代表被排除分子都没有活性。下一步用已校准预算接实际 WEE1 共晶查询与 E024
精细三维匹配，测量其额外耗时与压缩；有已知活性／非活性标签后再做富集评估。
本次源名称碰撞的立体化学 QC 仍独立保留，不因索引验收通过而自动解决。

## 可选参数

```bash
# 自定义输出和 CPU 数量（输出目录必须是全新的）
E033_THREADS=20 E033_OUTPUT=/mnt/local/hand/yuzhang/aidd/e033-review-2 \
  bash scripts/run_e033_library_acceptance.sh --query-count 16

# 已有独立查询描述符时替换库内查询
bash scripts/run_e033_library_acceptance.sh --queries /absolute/path/queries.npz
```

外部 NPZ 必須包含原始、未经 z-score 的 `vectors`（n×60）和 `names`（n）；可选
`exclude_molecule_ids`（n）指定需要排除的同分子。必须按现有相同 USRCAT 生成方法
制备，并提供来源。禁止用不明模型产生的任意 60 维向量替代。

`--metadata-only` 可以跳过大文件哈希，但报告只标记元数据通过，不能当作完整
字节验收。首次正式运行建议使用默认完整检查。精确参考扫描在两种模式都执行。

## 方法来源

采用 FAISS 官方文档中的搜索参数／速度准确度权衡思路；参数选择和阈值属于
本项目 E033 协议，不是官方保证：
https://github.com/facebookresearch/faiss/wiki/How-to-make-Faiss-run-faster
https://github.com/facebookresearch/faiss/wiki/FAQ
