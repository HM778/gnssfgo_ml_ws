继续分析工程结构和数据流。

工程结构清楚了。现在读取关键文件：gnssfgo 侧的残差计算与桥接输出，以及 Transformer 侧的分析器。

现在深入读取 gnssfgo 侧的核心文件。

桥接层结构清楚了，注意到几个可疑点（缺省值填 0）。继续读 gnssfgo 主节点代码。

主节点代码读完，已经发现几个可疑点（`operator[]` 缺省插入、iono 参数置零等）。继续读因子图模型和伪距/载波相位因子定义。

伪距/载波因子已读。现在读因子图模型核心 `SelfAdjustedFactorGraph.hpp`。

因子图模型读完，发现了 `CPresidualsUpdate()` 的量纲问题。现在转向 OSQA/Transformer 侧，读分析器核心代码。

融合逻辑清楚了：`weighted` 模式下 suspect 意味着加权分落在 [0.2, 0.6)。继续读三个分析器。

Transformer 与图分析器都有系统性压分嫌疑。继续读时序分析器和桥接/读取代码。

反馈链路里发现一个关键问题（`_prn_to_sat_id` 的卫星编号对不上）。继续核对 gnssfgo 伪距因子构建与校正量计算，确认导出残差与实际因子是否一致。

核对卫星编号规则和 `getQualityScore` 实现。

卫星编号规则确认了：GLONASS/Galileo/BeiDou 的 `sat_id` 回传映射是错的。再看看有没有实际运行的输入/输出数据可以佐证。

没有历史运行数据。我用一个小的数值仿真验证 Transformer/图分析器的系统性压分假设。

Transformer 部分已证实：健康卫星 q_t ≈ 0.275（熵项归一化问题）。修复仿真脚本里的一个广播 bug 继续验证图分析器。

分析完成。先给结论，再列完整的原因清单和验证方法。

# 结论（TLDR）

**全部卫星被判 suspect 的主因在 Transformer/OSQA 侧的评分设计，而不是 gnssfgo 的伪距/载波残差计算方式。** 我用一个“完全健康”的合成历元（SNR 35–50dB、仰角 15–80°、伪距残差 1.5m、锁定稳定）按 OSQA 的实际代码逻辑做了数值仿真，复现了你观察到的现象：**15 颗健康卫星中 14 颗被判 suspect**（q_final 均值 0.486，而 trusted 门槛是 0.6）。即使已经加了 `--permissive`，结论不变。

核心机制：**健康的卫星彼此特征相似 → 注意力接近均匀分布 → 熵接近最大值 → `score_entropy = 1 - H/log(N-1) ≈ 0.04`**，把每颗卫星的 Transformer 分数压到 ≈0.27。也就是说，这个指标惩罚的恰恰是“正常”状态，所有卫星无差别饱和在 suspect 区间，失去区分度。

---

# 一、主因：OSQA/Transformer 侧（评分设计问题）

### 1. Transformer 熵项结构性偏低（影响最大，几乎是决定性的）

`transformer_analyzer.py:310-318`：

```
score_attention = clip(attended_by × N/2, 0, 1)   ≈ 0.5（均匀注意力时，正常）
score_entropy   = 1 - H/log(N-1)                  ≈ 0.04（健康卫星必然接近均匀注意力）
score_memory    = clip(记忆库相似度)               （初期记忆库为空 = 1.0，无贡献）
quality = (三者几何平均)^(1/3)                     ≈ (0.5×0.04×1)^(1/3) ≈ 0.27
```

仿真结果：15 颗健康卫星的 q_transformer 全部落在 0.27 左右（entropy_norm ≈ 0.96）。**随机投影的注意力 logits 标准差只有 ~0.3（温度 1.3 时），对 N≈15 颗卫星做 softmax 本来就近乎均匀**，所以熵项对所有卫星都一样低——好的和坏的卫星都拿 ~0.27。

**可立即验证**：看 `osqa_output.jsonl` 里每颗卫星的 flags，如果都带 `high_entropy`（`transformer_analyzer.py:323-324`，熵 >0.7 触发），即证实此条。

### 2. 图分析器的一致性误差与特征尺度不匹配

`graph_analyzer.py:340`：`q_graph = exp(-consistency_error / T)`，T=0.9（permissive）。`consistency_error` 是原始特征与两层 GCN 聚合特征的 RMS 差（`graph_analyzer.py:327-328`）。z-score 归一化后，卫星间特征天然相差 ~1 个标准差量级，GCN 聚合（self_weight=0.5）必然把自己拉向邻居 → RMS 误差 ~0.5–1.2 → q_graph ≈ 0.3–0.6，系统性中等偏低（我仿真里 q_g ≈ 0.57）。误差尺度（特征空间）和温度（0.9）根本不在一个标定体系里。

另外 `graph_analyzer.analyze()` 的 `pseudorange_residuals` 参数**传进去了但从未使用**（`graph_analyzer.py:271`，整个函数体没引用）。

### 3. trusted 门槛在当前融合结构下几乎不可达

`quality_fusion.py:284-289`：trusted 需要 q_final ≥ 0.6；weighted 融合（`quality_fusion.py:252-254`）下 `q_final = 0.45·q_t + 0.40·q_g + 0.15·q_tmp`。当 q_t 被钉死在 0.275 时，即使 q_g=q_tmp=1.0，q_final 最高也只有 0.674；实际 q_g≈0.57、q_tmp≈0.9 → q_final≈0.49 → 全体 suspect。**0.6 这个门槛是在 q_t 正常分布假设下调的，与熵项 bug 叠加后变成“永不 trusted”。**

### 4. “在线学习”实际上是自锁死循环

记忆库只收 `quality_final ≥ 0.7` 的好样本（`run_analyzer.py:396`），而 q_final 永远到不了 0.7 → `good_samples` 永远为空 → `memory_similarity` 恒为默认 1.0（`transformer_analyzer.py:280`）→ 记忆增强通道完全没参与。这不是压分原因，但意味着 Transformer“学不到东西”，标注质量也无法随时间改善。

### 5. 两个死特征

- `run_analyzer.py:192` 调 `extract_batch(observations)` 时**没传 `satellite_histories`** → 特征[7] 仰角变化率恒为 0（`feature_extractor.py:178-187`）。
- 特征[5] 载波残差恒为 0（见下文第三节第 1 条，gnssfgo 键名不匹配导致 OSQA 永远收到 0）。
- 8 维特征实际只有 6 维在动，且均含较大公共分量，进一步削弱区分度。

---

# 二、gnssfgo 侧：导出数据的问题（不是全部 suspect 的主因，但有几处实际错误）

### 1. 载波相位残差通道完全失效（键名不匹配）

- `SelfAdjustedFactorGraph.hpp:626-680` `CPresidualsUpdate()` 只往 `residuals[sat]` 里写了 `"dop_cp"` 一个键；
- `SelfAdjustedNode.cpp:582-595` 导出时 `res = factor_graph.residuals[sat]`（只有 dop_cp）+ `res["psr"]`；
- 而 `transformer_bridge.hpp:429-432` 和 OSQA 侧读的键是 **`tr_dd_pr` / `tr_dd_cp`** —— 不存在，恒导出 0.0。

结果：OSQA 的 `carrier_residual`（`run_analyzer.py:337`）永远是 0，载波残差特征和 `quality_fusion.py:258-268` 的 dop_cp 硬惩罚通道从未生效。**你精心做的 TR-DD-CP 残差反馈，OSQA 一位都没收到。**

### 2. dop_cp 残差本身量纲是错的

`SelfAdjustedFactorGraph.hpp:674-677`：`delta_cp_obs` 单位是**周**（cp 是周），`doppler_move = Δt × |v_receiver|` 单位是**米**（还是接收机速度的模长，没有投影到 LOS，也没有扣除卫星运动），`dop_cp_residual = 周 − 米`，物理意义不成立。另外 `sat_cp_const_l1[sat]` 用 `operator[]` 访问，未初始化的卫星会插入 0.0，使 delta_cp_obs 直接变成整段载波变化量（数百周）。

顺带：`SelfAdjustedFactorGraph.hpp:646` `obs->cp[l1_idx]` 没检查 `l1_idx >= 0`，单频卫星会读 `cp[-1]`，是未定义行为——**这很可能与你 todo.md 里记的“100 秒左右崩溃”有关**，建议单独排查。

### 3. `psr_residual_l1` 实际上是 L2 的残差

`buildCorrectedPseudorangeMeasurements` 对每颗卫星依次 push L1、L2 两条测量（`gnss_comm_extra.h:196-239`），而 `SelfAdjustedNode.cpp:549-578` 用 `psr_residual_map[sat] = residual` 逐条覆盖，map 里最终留下的是 **L2 的后验残差**，却以 `_l1` 的名字导出。伪距后验残差本身与 `pseudorangeFactor` 的模型一致（sagnac/钟差/TGD/对流层符号都对，电离层参数两边同样置零），方向没错，但标签错了。

### 4. sv_info 缺失时导出全 0

`transformer_bridge.hpp:396-415`：卫星不在 `sv_info_map` 或 `avaliable=false` 时，elevation/azimuth/sat_pos 全填 0。启动初期 ENU 参考点未建立前（`SelfAdjustedNode.cpp:353` 才写入 sv_info），**所有卫星都以仰角 0°、方位角 0°、卫星位置 (0,0,0) 导出**。OSQA 没有任何“仰角=0 无效”的判别，这些样本会污染前若干个历元的特征统计和记忆库。

### 5. 次要问题

- `SelfAdjustedNode.cpp:156` 读回传用的是**不带 time_frame 匹配**的 `readLatestQualityScores(path)` 重载（`transformer_bridge.hpp:518-531`），直接取文件最后一行——OSQA 处理滞后时，分数会套到错误的历元上（桥接头文件里专门写的 time_frame 匹配版本没被用）。
- 导出的后验残差未乘因子置信度（因子内是 `(est-obs)×1/σ²`），OSQA 拿到的是裸米值，两者不一致但影响温和。
- `receiver_enu` 导出的实际是 SPP 的 `latest_pos_enu` 而非注释所说的“优化后”。

---

# 三、回传链路：卫星编号映射错误（会误伤，但不是全部 suspect 的原因）

`gnssfgo_data_reader.py:183-207` `_prn_to_sat_id()` 对 R/E/C 卫星**只返回系统内 PRN 号**（R03→3、E01→1、C14→14），而 gnss_comm 的实际编号是分段偏移的（`gnss_utility.cpp:48-66`）：GPS 1–32、**GLONASS 33–59、Galileo 60–97、BeiDou 98–160**。

后果：
- 输出 JSON 的 `satellites` 以 `str(sat_id)` 为键，G03 和 R03 都映射成 `"3"`，**互相覆盖，非 GPS 卫星大量丢失**；
- gnssfgo 拿到的 sat_id 与 `obs->sat` 对不上，质量分数会**套到错误的卫星上**（R03 的分数加到 GPS G03 头上），非 GPS 卫星大多落到默认 1.0。

这解释不了“全部 suspect”（那是 OSQA 判的），但等 OSQA 修好后，这条会让权重反馈完全错位，必须一起修。

---

# 四、不改代码的验证路径

1. **看 `osqa_output.jsonl`**：如果每颗卫星 `details.transformer` 都在 0.25–0.35、flags 都含 `high_entropy`，即坐实主因是 Transformer 熵项，与 gnssfgo 残差无关。
2. **对比三个分量**：`graph` 是否也在 0.3–0.6 徘徊（图分析器尺度问题）；`temporal` 是否接近 1.0（时序正常，说明 gnssfgo 导出的观测序列本身是平稳健康的——反过来说明输入数据没问题）。
3. **OSQA 终端加 `--debug`**，看每历元 `Avg quality` 与各分析器均值（`get_long_term_stats()`）。
4. 用 `test_osqa.py` 离线回放一段已知的良好数据，确认同样复现全 suspect——与 gnssfgo 运行状态完全解耦。

# 五、修复方向建议（供后续参考，未改任何代码）

1. **熵项重标定**（收益最大）：把 `1 - H/log(N-1)` 换成相对基准的度量，例如以均匀注意力熵为“满分基线”做 `clip((H_max - H)/span)` 的反向映射，或干脆只用 `attended_by × N/2` 与记忆相似度，去掉熵项。
2. **图分析器校准**：`consistency_error` 除以邻居特征差异的自然尺度（如邻居间 RMS 距离），或把温度调到与实测误差同量级。
3. **trusted 门槛**与融合结构一起重调（q_t 修复后 0.6 或许可达，但建议先量化 q_final 分布再定）。
4. **对齐载波残差键名**（`dop_cp`→`tr_dd_cp` 或反向），并修正 dop_cp 的量纲（乘波长转米、LOS 投影、扣卫星运动），修掉 `cp[l1_idx]` 越界。
5. **`_prn_to_sat_id` 直接复用 gnssfgo 导入时写的正确 `sat_id`**（输入 JSONL 里本来就有，输出时别重算）。
6. 让 gnssfgo 读取时使用带 time_frame 匹配的重载，消除滞后错位。

一句话回答你的二选一问题：**两边都有问题，但“全 suspect”的直接原因是 Transformer 的评分设计（熵项 + 图一致性尺度 + 0.6 门槛的叠加），gnssfgo 的残差计算本身方向正确，主要缺陷是载波残差通道断链（恒 0）和若干标签/编号错误。**

---

# 六、2026-08-29 追加：天顶卫星 unreliable 的诊断（基于 16:14 的 osqa_input/output.jsonl，81 历元）

### 结论
天顶高仰角卫星（E88 73.7°、E92 67.6°、R53 75.4°、C28 71.6°、G20 78.1°）全程 unreliable 的直接原因
**不是**优化位置异常，而是：

1. **该数据由旧二进制产出**（进程内存中仍是旧 dop_cp 公式，磁盘 15:59 的新二进制未被加载）。
   dop_cp 残差呈"每星恒定偏置"（E88≈+237、E92≈−191、R53≈+445、C19≈+746、E78≈−694、C28≈+496，
   单位是周−米混合，不可解释），且恰好在该星载波相位对首次可用时出现。
2. OSQA 按米做硬惩罚（0.5m 起罚、9.5m 满罚 0.9）→ 所有 dop_cp≠0 的卫星 q×0.1 → q≈0.06~0.095 unreliable。
   天顶卫星载波跟踪最连续、prev-cp 对最早可用 → 最先被压死且永不恢复。
3. 铁证：dop_cp=0 的卫星全部 trusted（含 7.6° 低仰角的 R40），dop_cp≠0 的全部 unreliable（含 73.7° 的 E88），
   与仰角完全无关；G20 在 ts=979→980 之间 lock_count 0→1、dop_cp 0→−135.4，同历元 q 从 0.69T 掉到 0.08U。
4. **处置**：重启 gnssfgo 加载新二进制即可。验证方法：新 osqa_input.jsonl 中 dop_cp_factor_residual |值| 应 < 1~2m。

### 顺带证实的两个真实问题
- **BeiDou PRN 字符串截断**：导出 prn 用 `sat % 100`，BeiDou sat_id ≥ 100 丢失百位
  （显示 "C19" 的实际是 sat 119 / BDS prn22）。OSQA `_prn_to_sat_id` 反推得到 116（−3），
  评分会套到错误的 BDS 卫星上；且 C28(→125) 与 E66(→125) 在输出 dict 的 key 上碰撞，
  C28 的评分整个丢失（本次数据 C28 输出条目为 0 即此原因）。
  **修复方向**：OSQA 直接复用输入里已导出的正确 `sat_id` 字段（输入 JSONL 中 C28→128 是对的），
  不要从 PRN 字符串反推；C++ 导出 prn 改用真实 PRN（sat − 系统偏移）。
- **FGO 优化位置确实偏离伪距几何（用户猜想成立但不影响可信度判定）**：
  用每历元伪距后验残差反演 r_i = c_sys + LOS_i·δ，得 |δ| 平均 8.4m（最大 17.6m）、
  相邻历元变化平均 3.3m（最大 11.7m），以 z/竖直分量为主。
  嫌疑集中在 TR DD CP 因子的模糊度处理（sat_cp_const 快照机制）把解拉离伪距一致解。
  这影响 FGO 定位精度本身；扣除公共偏差后天顶星 psr 残差仅 −2~−3m、低仰角星 +4~+9m（未建模电离层），
  不足以触发任何不可信判定（融合对 psr 无硬惩罚）。建议单独排查模糊度/锚定问题。

---

# 七、2026-08-29 追加 2：上述问题的修复与根因补充

### 关键新发现：观测流为合成数据, cp/doppler 与星历几何不自洽
- 各卫星 cp/doppler 速率彼此差异极大（G20 ≈1.9 m/s, C19 ≈439 m/s）, 不符合 MEO 轨道动力学;
- 但每颗卫星自身满足 **Δcp = −doppler·Δt（比值 1.0000）**, 即 cp↔doppler 自洽;
- 数据中存在真实毛刺（E66/E73 最大 38~42 周的相干性跳变）。
- 结论: 一切依赖星历几何的 cp 检验（含旧 sat_cp_const 机制与几何预测法）在该数据上都不可用,
  只有不依赖几何的 cp↔doppler 相干性（TDEC）检验才同时适用于合成与实测数据。

### 已实施的修复
1. **cycleSlipDetect → TDEC 检验**（SelfAdjustedFactorGraph.hpp）:
   预测改为 Δcp_pred = −doppler·Δt（相邻历元多普勒均值, 梯形积分）, 阈值 2 周;
   彻底移除 sat_cp_const_l1/l2 / sat_const_residual_l1/l2 快照机制
   （旧机制快照失稳 + 预测漏卫星运动 → 每 3~5 历元误报周跳、锁定反复清零,
     这就是 lock_count 序列 1,2,3,0,1 的来源, 也是 CP 因子反复失锁的根源之一）。
2. **CPresidualsUpdate → TDEC 相干性残差**: dop_cp = λ×(Δcp + doppler_mean·Δt),
   健康链路厘米级, 单位米, 与 OSQA 硬惩罚量纲一致。
   注: 本数据 cp 约定为 Δcp = −doppler·Δt（相应用 +ρ/λ 符号）, 与 TRDDCP 因子的
   DD 符号约定一致, 已实测验证。
3. **BeiDou PRN 截断修复**（transformer_bridge.hpp）: PRN 字符串改用 satsys 反解的
   系统内真实 PRN, 不再用 sat % 100。
4. **OSQA sat_id 直通**: RawObservation.sat_id → fusion(SatelliteQuality.sat_id) →
   输出 JSON 的 satellites key 与 sat_id 字段, 全程不再从 PRN 字符串反推
   （_prn_to_sat_id 仅保留为 CSV 回放回退）。

### 修复后重放验证（81 历元真实数据, TDEC 重算 dop_cp）
- dop_cp 分布: 中位数 12.4cm, p90 33.3cm, max 9.66m（毛刺被正确标出）;
- 天顶卫星 E88/E92/R53/G20/C28 全程 trusted（0.73~1.00, 仅在数据毛刺历元瞬时 suspect 后恢复）,
  修复前它们被钉死在 0.05~0.10 unreliable;
- 全星座 trusted 19~30 颗 / unreliable ≈0（修复前 trusted 3~17 且持续劣化, unreliable 15~27）;
- sat_id 回传 33/33 一致、无碰撞, C28 评分不再丢失（修复前被 E66 覆盖）。

### 遗留事项
- 位置拉偏 |δ|≈8m 的剩余部分: 合成观测流与星历几何的固有不一致（DD 后残余 ~2.4m
  被 σ≥0.2m 的 CP DD 因子放大）+ CP 因子权重下限 max(0.5, ...) 在 OSQA 压权时削弱伪距锚定。
  实测数据上 CP 因子观测将回归自洽, 预计显著缓解; 若仍明显, 再排查锚定/模糊度。

---

# 八、2026-08-29 追加 3：伪距/载波因子权重量纲审计与优化策略重平衡

### 权重审计结论（修复前）
| 因子 | 归一化 | 有效 σ | 每米乘子 | 数量 |
|---|---|---|---|---|
| 伪距单差 | min(1, 1/σ_eff²), σ_eff=√(psr_std/0.16)/q | ≈1.4m | ≈0.5 | 120 (24/历元×5) |
| 双差伪距 | conf/√(Σpsr_std²) | ≈0.64m | ≈1.6 | 与TR对同数 |
| TR DD CP | /max(0.20, avg(σcp)) × tr_score | 0.2m(下限淹没) | 2.5~5 | 20 |
| 多普勒 | var 作为乘性置信度 | — | 0.5~1 | ~4 |

- getVarofCp_ele_SNR 因 useEleVar 未定义走了 0.0001·√(1/var) 分支（毫米级物理量），
  被 max(0.20,·) 下限完全淹没 → TDCP 精度被人为劣化 20 倍以上;
- tr_score 的 0.5 下限使低可信卫星仍保留一半 CP 权重（可信对比度不足）;
- 伪距 confidence 公式把权重当 σ 用（方向倒置），被 1.0 钳位掩盖;
- 多普勒 quality_scale 取倒数（低质量→约束增强），同样被钳位部分掩盖;
- CP=20 是运行时控制器被 opt_time≈56ms 压到下限（MaxTRFactorNum 收缩）。

### 已实施的重平衡
1. getVarofCp_ele_SNR 改为高程模型 σ（天顶 1.4cm / 15° 4cm, 单位米）;
2. TRDDCP 双差 σ 按独立观测平方传播（≈2× 单差）, 下限 0.10m
   （实测校准: 0.03 收敛段瞬态过大 z≈6m, 0.10 收敛段优于旧配置）;
3. tr_score 下限 0.5→0.1: 可信卫星 (q≈1) CP 全额权重, 低可信按比例降权;
4. 多普勒 quality_scale 方向修正 (×avg_quality);
5. launch: max_psr_factors_per_epoch 24→12, min_tr_factor_num 20→60
   → CP 20→60, 伪距 120→60, opt_time 54ms 在 50ms 预算附近。

### A/B/C 无头重放实测（同一 bag 前 120s, 二阶差分抖动）
| 配置 | 收敛段\|d2\| | 稳态\|d2\| | 稳态 x/y/z std | CP数 |
|---|---|---|---|---|
| A 旧 (σ0.20/CP20/psr120) | 2.377 | 1.440 | 1.04/0.28/1.85 | 20 |
| B 强CP (σ0.03) | 3.251 | 0.762 | 0.32/0.32/0.76 | 60~129 |
| **C 采用 (σ0.10)** | **1.733** | **0.755** | **0.32/0.32/0.75** | **60** |

结论: 稳态平滑度提升约 2 倍, 收敛段也优于旧配置; TDCP 高精度约束进入稳态后
成为相对轨迹的主导约束, 伪距 (数量减半) 负责绝对锚定, 与用户预期一致。

---

# 九、2026-08-29 追加 4：遮挡段退化感知修复 + 反馈链路 time_frame 断链修复

### 现象与根因
最后 90s (300-394s) 遮挡段: SPP |d2| 抖动 1.0→3.5m, p90 伪距残差 5-10→19.2m,
但 OSQA mean_q 反而升到 0.824 —— 三个分析器全为相对比较机制
(与同伴比/与自己 EMA 比/与邻居比), 共性渐变退化随归一化统计一起漂移, 集体失明;
而绝对证据 (伪距后验残差/SNR) 只存元数据不参与打分。

### 修复
1. **融合层新增绝对健康度惩罚** (quality_fusion.py, 不随群体/历史自适应):
   pen_res = clip((|psr残差|−3)/9, 0, 0.9); pen_snr = clip((30−SNR)/12, 0, 0.9);
   q_final ×= (1−pen_res)(1−pen_snr); 新标记 high_psr_residual / low_snr。
2. **反馈链路 time_frame 断链修复** (gnssfgo_data_reader.py):
   write_quality_result 此前不写 time_frame 字段, 而 gnssfgo 端用
   readLatestQualityScores(path, time_frame) 按历元匹配 → 匹配永远失败 →
   getQualityScore 恒为默认 1.0, **质量反馈在优化器侧从未真正生效**
   (改用匹配重载时引入的回归)。现输出补写 time_frame = timestamp × 10。

### 全 bag 在线验证 (gnssfgo + OSQA 同跑, 386 历元)
- OSQA 评分: 退化段 mean_q 0.824→0.616 (最低 0.391); 逐星区分清晰:
  受遮挡星 E-1/G07/R-24 → 0.03-0.15 (high_psr_residual/low_snr),
  健康星保持 0.75-0.89; 健康段整体仅下移 ~0.08, 无新增误报。
- **公里级单历元发散**: 无OSQA基线 8 处 (144/145/199/200/313/314/317/318),
  反馈修复后仅剩健康段 2 处 (144/145), **遮挡段发散全部消除**。
- 未竟事项: t≈144-145 健康段的发散在两组运行中都出现, 与遮挡无关,
  疑似窗口锚定/个别 CP 对初始化问题, 需单独排查;
  src/c099.csv 为追加式输出, 已混入多次验证重放的行
  (本次验证行 trddcp_factor_count=60, 可据此过滤), 需要时可整文件轮转清空。

---

# 十、2026-08-29 追加 5：输出跳跃与频率不均的修复

### 跳跃根因
t=144 处接收机钟跳步(stepping)使所有卫星 Δcp 同步偏移 → TDEC 检测集体误报周跳
→ 锁定计数集体清零 → TDCP 约束瞬间真空 → 优化器收敛到公里级错误盆地
(opt_time 飙至 300-440ms, 协方差却不大), 下一历元自愈。

### 修复 (SelfAdjustedFactorGraph.hpp / SelfAdjustedNode.cpp / base_factorgraph.hpp)
1. TDEC 周跳检测加共模抑制: 扣除本历元全体卫星 TDEC 残差的中位数
   (按历元+频段缓存), 钟跳步属接收机行为不再误报为周跳;
2. 求解后发散防护 guardLatestStateAgainstDivergence(): 最新状态与 SPP 初值
   偏差 >50m 或非有限时回退为 SPP 初值并 WARN (实测捕获一次 2318m 跳变);
3. max_num_iterations 50→20: 滑窗热启动通常数次迭代即收敛,
   opt_time max 从 300-440ms 降到 171ms, 消除挤占下一历元处理窗口的尖峰。

### 验证 (全 bag + OSQA 在线, 394 历元)
- 输出间隔 393/393 全部 1s (无丢历元/重复);
- |d1|>50m 跳变 8 处 → 0 处 (防护捕获 1 次);
- 稳态 |d2| 中位 0.865m (持平), 遮挡段 3.41m (与 SPP 的 2.95m 接近,
  含真实车辆动态; FGO 的价值在于消除公里级发散而非降低该段噪声)。
- 遗留: t=145 仍有部分锁定清零 (15/29, 共模中位数未完全吸收钟瞬态);
  t≈144-145 的发散被防护兜住但根因未完全消除; 用户满负载环境 (rviz+--vis)
  若仍偶发丢历元, 可减小 max_psr_factors_per_epoch 或关闭 OSQA 可视化。
