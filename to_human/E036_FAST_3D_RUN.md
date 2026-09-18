# E036：并行 chunk 与全库 3D 搜索

E035 百万验证已通过（用户提供的报告，未在本机复核远端原始文件）：
1,000,000 构象 / 952,402 来源分组分子 / 6,000,000 对比；误差 0、分配差异 0，
构象及分子排名一致。QT9 E031 1.4476 秒（历史 refine 的 4.63%），824 1.9092 秒
（2.73%）。E031 保持注释用途。已导出 39 / 41 个姿态，仍待人工检查。

## 优先执行：实际大库搜索

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e036_fast_3d.sh
```

默认输出：`/mnt/local/hand/yuzhang/aidd/e036-fast-3d/run-w16-c500-r64-v1`。
默认原始 E034：`/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/recheck-20260918-163847-316118`。
运行两条已校准的 WEE1 查询，复用现有全库 FAISS 索引，nprobe 128，预算 10,000。
精确重算候选的描述符距离，随后 Gaussian coarse / refine 和 batched E031。
这不是任意新靶点已验证的搜索服务；本阶段验证全库搜索执行链与已保存结果等价。

调度为 16 workers，coarse chunk 500，refine chunk 64；原 E034 为 2000 / 250。
评分方程、5000 Top-N、512 seed cap、候选预算均保持不变。原 coarse 只有5个任务，
本配置有20个；refine 约99–107个任务，有望降低尾部等待，需要工作站实测。
任务队列最多保留 2 × workers 个 future。默认查询间顺序执行，索引只加载一次。

每条查询完成后，比较所有 Gaussian 数组（包括变换矩阵），以及 E031 分数、分配
和排名；不一致则报错。新的检索候选也必须与 E034 保持一致。完整性检查、索引加载、
最终等价性比较单独处理；不运行全库精确参考扫描或百万压力测试。
`execution_seconds` 包括候选一致性检查和输出写入，不能称为服务 SLA。

发送新目录中的 `report.md` 和 `report.json`。观察总查询耗时和 coarse/refine wall time，
不只观察单块秒数。并行加速比目前未在真实工作站上测量。

中断后可恢复：

```bash
bash scripts/run_e036_fast_3d.sh --resume
```

完整查询复用会保留原始计时；部分 chunk 复用的查询标记 `fresh_compute: false`，
不能用来判断全新查询的耗时。需要纯计时时指定新的 `E036_OUTPUT`。

## 并行大规模验证（不必重跑已通过的百万验证）

E035 runner 新增 `E035_WORKERS`（默认8）和 `E035_CHUNK_SIZE`（默认512），
使用 spawn 进程、各自 mmap reader、稳定顺序合并、原有哈希断点。
每10块显示进度和 ETA，`scale/progress.json` 每块更新。ETA 不包含最终全局排名检查。
worker 累计 scoring seconds 与实际 wall seconds 分开报告。

如需衡量并行提速，可在空闲工作站上运行固定100k样本的串行/8进程对照：

```bash
E035_WORKERS=1 E035_SCALE=100000 \
E035_OUTPUT=/mnt/local/hand/yuzhang/aidd/e035-review-and-scale/scheduler-serial-100k-v2 \
bash scripts/run_e035_review_and_scale.sh

E035_WORKERS=8 E035_SCALE=100000 \
E035_OUTPUT=/mnt/local/hand/yuzhang/aidd/e035-review-and-scale/scheduler-parallel-100k-v2 \
bash scripts/run_e035_review_and_scale.sh
```

比较 `scale.execution.scale_wall_seconds_this_invocation`，确认两次均无复用，
样本哈希和所有等价性门槛相同。运行顺序有缓存影响，单组对照属探索性测量。
8个进程是保守起点，不承诺8倍加速；内存、磁盘、分子复杂度影响吞吐。
旧 `run-1m-v1` 保留；代码及调度改变后不要拿旧协议目录继续新实验。

## 下一步判断

如更细分块收益小，需基于 Gaussian worker 的计时分析 seed 搜索、距离矩阵和对象读取，
再进行评分等价的内核优化。减少候选或 seed cap 属另一类实验，必须另测覆盖损失。
分子级 Top100 交集仅6 / 9，不证明 E031 优于 Gaussian；两种评分差异姿态仍待检查。
