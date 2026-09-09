---
name: setup-data
description: Download the TiPhD benchmark test dataset (~hundreds of MB) from ModelScope (fansailing/TiPhD_test_data) to <cwd>/data/. Idempotent via a marker file (.tiphd_data_ready). Use this skill when the user says "download data", "setup data", "下载数据", "TiPhD data", or before invoking any algorithm skill for the first time.
---

# setup-data Skill

> 从 ModelScope 下载 TiPhD 测试数据集到本地。**幂等**——已下载会跳过。
> 默认下载到 `<cwd>/data/`；可通过环境变量 `TIPHD_DATA_DIR` 覆盖。

---

## 一、本 skill 做什么

收到用户指令后，自动完成：

1. **幂等检查**——查 `<target>/.tiphd_data_ready` 是否存在
2. **下载**——从 ModelScope `fansailing/TiPhD_test_data` snapshot 到本地
3. **完整性验证**——查 25 个关键文件是否齐全（5 癌种 scRNA + bulk + gold standard）
4. **写 marker**——下次 setup-data 看到 marker 立刻跳过

---

## 二、关键路径

- 下载脚本：`.claude/skills/_toolkit/download_data.py`
- 默认目标：`<cwd>/data/`（可用 `TIPHD_DATA_DIR` 覆盖）
- 数据源：https://www.modelscope.cn/datasets/fansailing/TiPhD_test_data/files
- 完整性 marker：`<target>/.tiphd_data_ready`

---

## 三、标准作业流程

### 步骤 1 — 幂等检查

```bash
SKILL_DIR="<path to .claude/skills>"
TARGET="${TIPHD_DATA_DIR:-$(pwd)/data}"
MARKER="$TARGET/.tiphd_data_ready"

if [ -f "$MARKER" ]; then
  echo "OK data already at $TARGET (marker exists), skipping"
  echo "   delete $MARKER to force re-download"
  exit 0
fi
```

### 步骤 2 — 下载

```bash
pip install modelscope   # 一次性
python "$SKILL_DIR/_toolkit/download_data.py" --target-dir "$TARGET"
```

实现：
```python
from modelscope import snapshot_download
snapshot_download(
    repo_id="fansailing/TiPhD_test_data",
    repo_type="dataset",
    local_dir=TARGET,
    revision="master",
)
```

首次可能下载几百 MB；ModelScope SDK 自动断点续传。

### 步骤 3 — 完整性验证

下载脚本自动检查 25 个关键文件存在：
- 5 癌种 scRNA（.rds 给 R，.h5ad 给 Python）
- AML 3 bulk + CRC 5 bulk + HCC 1 bulk + GC 1 bulk
- 5 gold standard CSV（含 `gold_standard_all_LUDA.csv` 这个**拼写错误**的文件名）
- 详见 `download_data.py` 的 `expected` 列表

### 步骤 4 — 写 marker + 报告

```bash
du -sh "$TARGET"
ls "$TARGET"
```

---

## 四、用户的可选参数

| 用户说法 | skill 行为 |
|---------|-----------|
| "setup data" / "下载数据" | 默认 <cwd>/data；完整 4 步 |
| "下载到 /scratch/data" | `TIPHD_DATA_DIR=/scratch/data python download_data.py --target-dir /scratch/data` |
| "force re-download" / "重下数据" | 删 `<target>/.tiphd_data_ready` 后重跑 |
| "verify data" / "验证数据完整性" | 跑 download_data.py，看是否走 "already ready" 分支 |
| "data 在哪" / "data path" | 报 `os.environ["TIPHD_DATA_DIR"]` 实际值；fallback 是 `<cwd>/data` |

---

## 五、注意事项

1. **数据放在哪？** 默认是用户**当前工作目录**下的 `data/`。**用户运行 skill 时所在的目录**就是下载目标。建议在项目根目录下运行（如 `cd ~/myproject && claude`）。
2. **LUDA 文件名拼写错误**：原项目文件名是 `gold_standard_all_LUDA.csv`（不是 `LUAD`）。skill 已适配，**不要改名**。
3. **scRNA 文件扩展名**：R 算法用 `.rds`，Python 算法用 `.h5ad`；两个都有下载。
4. **下载时间**：取决于网络，通常 1–10 分钟。
5. **磁盘**：5 癌种所有数据几百 MB 到 1GB 级别。
6. **重下**：删 `<target>/.tiphd_data_ready` 即可，modelscope 会断点续传。
7. **数据源迁移**：本项目历史数据曾在 HuggingFace（`SailingFan/TiPhD_test_data`），现已迁到 ModelScope（`fansailing/TiPhD_test_data`）。**只用 ModelScope**。

---

## 六、失败兜底

| 现象 | 处理 |
|------|------|
| `ModuleNotFoundError: modelscope` | `pip install modelscope` |
| 网络超时 | 重试；国内用户可设 `MODELSCOPE_CACHE` 到本地路径加速 |
| 部分文件缺失 | 删 marker 重跑；modelscope 会重新拉 |
| 磁盘满 | 换目标目录（`--target-dir`）；先清旧数据 |
| 权限拒绝 | 换 `TIPHD_DATA_DIR` 到有写权限的目录 |

---

## 七、不在 skill 范围内

- 上传数据到 ModelScope（不在 skill 职责内）。
- 修改数据内容（只下载，不改）。
- 装 modelscope（让 setup-env 处理）。
