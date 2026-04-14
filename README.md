# 3状態 Bass 拡張モデルによる日本 COVID-19 新規陽性者数のモデリング

マーケティング動態の **3状態 Bass 拡張モデル**（SIRS 同型）を日本国内 COVID-19 新規陽性者数（2020-2023/5、週次）に直接適用し、複数波をベイズ推定で同時分離する試み。

## 背景

3状態 Bass 拡張モデル（加入・離脱・再加入）は疫学の **SIRS モデルと数学的に同型**：

| Bass | 疫学 | COVID |
|---|---|---|
| s(t) | Susceptible | 未感染者 |
| a(t) | Infectious | 感染中 |
| x(t) | Recovered | 既感染・回復 |
| α | 再感染率 | 免疫減衰後の再感染 |

## モデル

$$
\begin{aligned}
\frac{ds}{dt} &= -(p_1(t)+q_1 a) s \\
\frac{da}{dt} &= (p_1(t)+q_1 a) s + \alpha (p_1(t)+q_1 a) x - \gamma_{rec} a \\
\frac{dx}{dt} &= \gamma_{rec} a - \alpha (p_1(t)+q_1 a) x
\end{aligned}
$$

時変 $p_1(t)$ は4波（Delta/BA.1/BA.5/XBB）のガウス摂動重ね合わせ：

$$p_1(t) = p_{1,\text{base}} + \sum_{k=1}^{4} A_k \exp\!\left(-\tfrac{1}{2}\left(\tfrac{t-c_k}{w_k}\right)^2\right)$$

観測は新規感染フロー $\text{new}(t) = M \cdot (p_1(t)+q_1 a) \cdot s$ を LogNormal 尤度で。

## 主要結果

| パラメータ | 推定値 | 解釈 |
|---|---|---|
| M | 5240万人 | 感受性人口（総人口の約40%） |
| q₁ | 27.17/年 | 基本感染力（R₀ ≈ 1.045） |
| α | 0.16 | 再感染率 |
| σ_log | 1.29 | log 残差 std |

- 全パラメータ **rhat ≤ 1.003** で収束
- **BA.1・BA.5・XBB の3波をガウス摂動として明確に分離**
- Delta 波は q₁ の内因的伝播に吸収され消失
- R₀ ≈ 1.04 → 各波は **p₁(t) 外部ショックで駆動**（変異株・接触増）

## ファイル

| ファイル | 内容 |
|---|---|
| `COVID_Japan_3State_Bass.ipynb` | メインノートブック |
| `newly_confirmed_cases_daily.csv` | 厚労省オープンデータ（全国日次新規陽性者） |

## 実装環境

- Julia 1.10
- Turing.jl（NUTS, target accept 0.9, max_depth 10, 1000×4 chains）
- DifferentialEquations.jl（Tsit5, abstol=reltol=1e-6）

## 限界と今後

- 報告率 ρ(t) 未補正（2022年後半の検査縮小期を過小評価）
- 変異株ごとの q₁(t) 時変化未導入
- 二重尤度（抗体保有率調査との併用）で M・α の識別性改善余地あり
- 4状態 SEIRS 拡張で潜伏期 E を追加する方向性

## 結論

3状態 Bass 拡張モデルによる COVID-19 モデリングは技術的に成立。主要な寄与は **マーケティング動態と感染症動態の数理的同型性をベイズ推定で実証** した点。

## データ出典

- 厚生労働省 新型コロナウイルス感染症対策推進本部: https://covid19.mhlw.go.jp/
