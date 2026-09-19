# E038：Prompt → 蛋白数据 / AF3 / 已校准3D搜索

本次接口使用兼容 OpenAI 的 Chat Completions，台式机 AF3 已安装（用户确认）。
LLM 负责生成结构化任务；本地适配器读取真实蛋白序列、准备 AF3 输入并调用配置的
AF3 runner。AF3 本身接受结构化输入，不直接接受自由文本 prompt。

## 明天执行顺序

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main

# 无网络、无API Key、无GPU的接口测试
bash scripts/run_e038_prompt_smoke.sh

# 大库3D性能/等价性批量测试（若尚未运行）
bash scripts/run_e037_workstation_suite.sh
```

前一个脚本输出 `e038-prompt-smoke/<时间>-<PID>/report.json`，明确标记为离线模拟，
不代表真实LLM/蛋白API/AF3已经通过。后一个脚本仍按E037协议测量真实大库。

## 配置真实服务

在库外保存本地配置，例如 `/mnt/local/hand/yuzhang/aidd/config/`。
以 `configs/llm/openai-compatible.example.json` 为模板，填写你实际的 `base_url`
（到 `/v1`）、`model`；`api_key_env` 保留为 `AIDD_LLM_API_KEY`。
Key只在终端环境中设置，不写入JSON或Git，不发送到对话。
本地兼容服务可用 `http://127.0.0.1:<端口>/v1`，无鉴权时省略 `api_key_env`。
若服务不支持JSON mode，可将 `json_mode` 设为false，本地仍执行严格JSON校验。

```bash
export AIDD_LLM_PROFILE=/mnt/local/hand/yuzhang/aidd/config/llm.local.json
read -rsp 'LLM API Key: ' AIDD_LLM_API_KEY
echo
export AIDD_LLM_API_KEY
```

复制 `configs/prompt-runtime.example.json` 为库外 `prompt-runtime.local.json`，
将 `af3_profile` 指向已有AF3本地配置。配置格式见
`configs/models/alphafold3.example.json`：`python`、`runner`、`model_parameters`、
`databases`必须是实际存在的绝对路径，并满足项目原有的AF3 profile检查。
本工具不会安装AF3、下载权重或更改其环境。已安装的AF3 pipeline负责MSA/templates，
不会由LLM伪造这些输入。可用实际安装版本填写profile的 `version`。

```bash
export AIDD_RUNTIME_PROFILE=/mnt/local/hand/yuzhang/aidd/config/prompt-runtime.local.json
```

示例中的3D配置保留参考seed生成器（`bounded_pair_seeds:false`）。E037通过后，
根据结果选择分块及该选项；变更配置时创建新计划，不能混入旧协议目录。

## 先测真实LLM，再准备蛋白/AF3输入

```bash
# 只测试LLM计划，不请求蛋白API，不运行科学计算
bash scripts/run_e038_prompt.sh --plan-only \
  --prompt '查找人类WEE1的UniProt和分辨率不高于3埃的PDB候选，准备全长蛋白的AF3输入，暂不运行预测。'

# 执行数据读取和输入准备；输出FASTA、来源证据、PDB候选、AF3 JSON
bash scripts/run_e038_prompt.sh \
  --prompt '查找人类WEE1的UniProt和分辨率不高于3埃的PDB候选，准备全长蛋白的AF3输入，暂不运行预测。'
```

脚本自动建立/复用用户 `workstation` 的 `Prompt AIDD` Project并激活。
存储根默认 `/mnt/local/hand/yuzhang/aidd/prompt-workspace`，可用
`AIDD_PROMPT_STORAGE` 更改；脚本支持 `--username`、`--project-name`。
已有Project可直接用底层CLI指定真实 `--db --user --project`，它会检查所有权和激活状态。

输出打印 `plan.json`、context和执行报告位置。完整任务在Project的
`runs/PROMPT-<id>/` 下；LLM请求和计划也记录到现有AI审计表。
计划只允许固定动作，不能包含shell命令、任意路径、模型编造序列或检索预算。
模型云端上下文只包含用户输入和工具说明，不自动上传本地库、蛋白序列或实验数据。

## 真实AF3运行

```bash
bash scripts/run_e038_prompt.sh --allow-compute \
  --prompt '读取人类WEE1的UniProt全长序列，准备AF3输入并运行一次AF3预测，seed为1。'
```

也可明确指定UniProt accession、1-based闭区间残基范围或CCD配体，如：
“读取P30291全长序列，准备含ATP配体的AF3复合物输入，先不要运行。”
CCD身份会向RCSB核对。当前支持单条蛋白链加CCD配体，暂不支持prompt编排DNA/RNA、
多蛋白复合物、自定义CCD/共价键或自动域选择。需要更多信息时模型应返回澄清问题。
UniProt多个匹配不会自动选取；未知/非标准残基会停止AF3准备，不会擅自替换。

`--allow-compute` 是本次调用是否允许启动AF3/3D计算的执行开关。默认仍可完成数据读取
和输入准备；遇到计算步骤会保存 `blocked` 及原因。模型输出的运行请求不能绕过开关。
成功预测会保存 `prediction-manifest.json`、结构SHA、AF3置信度；它不是结合活性或
可靠pose的验收结论。失败尝试保留，新尝试写入独立目录。

## Prompt调用3D搜索

```bash
bash scripts/run_e038_prompt.sh --allow-compute \
  --prompt '使用已校准的WEE1 QT9查询对现有大库执行3D搜索和Gaussian精筛，E031只作注释。'
```

可选查询：`wee1_qt9`、`wee1_824`、`wee1_both`。当前是已校准WEE1链路的prompt入口，
不能把其他蛋白名自动替换为WEE1。新靶点仍需构建/验证配体query和受体结构。
LLM不改变10000候选预算、Top-N或seed cap，也不把E031升为筛选排名。
真实搜索保持与E034证据等价性检查。

## 恢复、回传

同一计划可用底层CLI恢复（将打印出来的真实ID/路径代入）：

```bash
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python -m aidd_agent.prompt_workflow run \
  --db /实际/registry/aidd.sqlite3 --user USR-实际ID --project PRJ-实际ID \
  --plan /实际/Project/runs/PROMPT-实际ID/plan.json \
  --runtime "$AIDD_RUNTIME_PROFILE" --allow-compute
```

恢复会校验计划、代码、本地配置和已完成产物；内容变更则要求新计划。
`ask`/上层脚本会创建新计划，不会自动发现旧计划。
报告区分 complete / blocked / failed / not_run；查看失败步骤目录的执行日志。
回传离线smoke报告、E037汇总，以及真实任务的 `plan.json` 和
`execution/report.json`（不要发送API Key、环境文件或敏感prompt）。

## 参考与验证边界

- [兼容接口依据：OpenAI Chat API](https://developers.openai.com/api/reference/resources/chat)
- [UniProt REST查询说明](https://www.uniprot.org/help/api_queries)
- [RCSB Data API](https://data.rcsb.org/)
- [AF3官方输入格式](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)

本实现使用AF3基本dialect v1蛋白/CCD格式以兼容已有安装。真实LLM、在线数据库、
GPU/AF3运行尚未在本机验证；离线和mock测试不替代明天的真实集成测试。
