# 1. はじめに　　

貯留関数法（角屋・永井, 1978）は，洪水流出解析で広く使われている。基礎式は，貯留量と流出量の関係を表す常微分方程式である。パラメータは通常，事前に決めた値に固定される。降雨中にオンラインで適応させる例は少ない。

数少ない例外が，ASFM（Kim et al., 1998, 2001）である。ASFMは，MMAE（複数モデル適応推定）を貯留関数法に取り入れた手法である。パラメータをリアルタイムで適応させる。実際の2流域のデータでも有効性が示されている。

ASFMの仕組みは，次のようなものである。

- 離散化したパラメータの候補を複数用意する（多重モデル）
- 各候補で計算した値と，観測値を比べる
- 当てはまりが良い候補ほど，ベイズの定理で重みを高くする。


本リポジトリではこれとは異なるアプローチを取る。対象は，星・山岡（1982）の分布遅延型の貯留関数モデルである。貯留量の水収支を，移動窓 $[t-\tau, t]$ 上の積分として書く。この形にすると，分布遅延をもつ同次系の一般安定性定理（Aleksandrov et al., 2023）がそのまま使える。この積分形への書き換えは，貯留関数でよく用いられる「その瞬間の流入・流出量だけで決まる形の微分方程式」の定式化とは異なる。

以上を踏まえ，本リポジトリの提案は次の３点である。

### (i) MRACによるオンラインパラメータ推定と安定性証明

MRAC（モデル規範形適応制御）により，本モデルは積分形の水収支で書かれている。上記の安定性の根拠を，分布遅延をもつ同次系の一般安定性理論（Aleksandrov et al., 2023）に結び付けることでLyapunov安定性証明を与える。未知パラメータをオンライン推定する。

### (ii) ポンプ場モデルとのリアルタイム結合

オンライン推定した流出量を，ポンプ場モデルへリアルタイムで結合する。ポンプ場は，可変速ポンプとポンプ井からなる。(Fecarotta et al., 2018から引用)

### (iii) sincカーネル遅延埋め込み + レプリケータ動力学によるピーク遅れ推定

降雨波形と推定流出量波形の間の，ピーク遅れ時間を推定する。sincカーネルによる遅延埋め込みと，レプリケータ動力学を用いる。


## 1.1 同次系としての定式化と安定性(引用)

出典: A. Aleksandrov, D. Efimov, E. Fridman, “Stability of homogeneous systems with distributed delay and time-varying perturbations,” *Automatica*, Vol. 153, 111058, 2023 の定理1(証明等は原論文参照)。

#### 同次性

$f:\mathbb R^n\to\mathbb R^n$ が*$r$-同次*(次数$\nu$)であるとは

$$
f(\Lambda_r(\lambda)x) = \lambda^\nu\,\Lambda_r(\lambda)\,f(x)
  \qquad(\forall x\in\mathbb R^n,\ \forall\lambda>0)
$$

<div align="right">式 (1)</div>
が成り立つことをいう。

$r=(r_1,\dots,r_n)$（$r_i>0$）は重みベクトル,

$$
\Lambda_r(\lambda):=\mathrm{diag}\{\lambda^{r_i}\}_{i=1}^n
$$

$$
r_{\min}:=\min_i r_i,\qquad \nu\ge -r_{\min}.
$$

#### 対象システムと仮定

$$
\dot x(t) = F(x(t)) + \int_{t-\tau}^{t} G(x(s))\,\mathrm{d}s
$$

<div align="right">式 (2)</div>

$x(t)\in\mathbb R^n$, $\tau=\text{const}>0$。$F,G$ は連続で, 重み $r$ に関して同一次数 $\nu$($\nu+r_{\min}>0$)の $r$-同次とする。

**仮定1**: 補助系(遅延なし対応系)

$$
\dot x(t) = F(x(t)) + \tau\,G(x(t))
$$

<div align="right">式 (3)</div>

のゼロ解が漸近安定(AS)。

#### 定理1

仮定1のもとで, (i) $\nu>0\Rightarrow$式(2)のゼロ解はAS, (ii) $\nu<0\Rightarrow$解は一様最終有界(UUB)。

#### なぜ成り立つか(証明の要点)

証明はLyapunov–Krasovskii汎関数(LKF)

$$
\tilde V(x_t) = V(x(t)) + \frac{\partial V(x(t))}{\partial x}^{\!\top}
    \int_{t-\tau}^{t}(s+\tau-t)\,G(x(s))\,\mathrm{d}s
    + \int_{t-\tau}^{t}\bigl[\alpha+\beta(s+\tau-t)\bigr]\|x(s)\|_r^{\mu+\nu}\,\mathrm{d}s
$$

<div align="right">式 (4)</div>

とおく。($V$ は仮定1が保証する次数 $\mu>2r_{\max}$ の同次Lyapunov関数, $\alpha,\beta>0$, $\|\cdot\|_r$ は同次ノルム)。Young/Hölderの不等式を用いて $\alpha+\beta\tau$ を十分小さく選ぶと, ある $h>0$ が存在して

$$
\dot{\tilde V} \le -h\,\tilde V(x_t)^{\,1+\nu/\mu}
$$

<div align="right">式 (5)</div>

が, ある領域上で成り立つことが示される。この領域, ひいては結論(i)(ii)の違いは指数 $1+\nu/\mu$ の大小関係から生じる。

- $\nu>0$ のとき $1+\nu/\mu>1$ であり, 式(5)は $\tilde V$ が小さい(原点近傍の)領域で成立する。指数が1を超える減衰は $\tilde V\to0$(したがって $x(t)\to0$)を導き, 局所的なASを与える。

- $\nu<0$ のとき $1+\nu/\mu<1$ であり, 同じ形の不等式は逆に $\tilde V$ が大きい領域で成立する。原点近傍では不等式が保証されないため, $x(t)\to0$ までは保証できない（UUB）。

不等式(5)の詳細な導出 (LKFの上下界評価, 微分の評価, 補題1による相殺)は原論文の定理1証明を参照。本モデルへの適用(仮定1の検証, 定理1/定理2の使い分け)は §1.2でまとめて行う。

#### 定理2(摂動系に対する一般安定性, 引用)

出典: 同上(Aleksandrov, Efimov, Fridman, 2023)の定理2。無摂動系(2)に, もう一段の分布遅延型積分項を加えた摂動系であり, 本モデル(§1.2, §1.3)で使用する。

$$
\dot x(t) = F(x(t)) + \int_{t-\tau}^{t} G(x(s))\,\mathrm{d}s
    + \int_{t-\tau}^{t} \mathcal R(s,x(s))\,\mathrm{d}s
$$

<div align="right">式 (5-1)</div>

を考える($F,G$ は式(2)と同一の仮定を満たす)。摂動 $\mathcal R:\mathbb R\times\mathbb R^n\to\mathbb R^n$ の大きさは, 各成分について

**仮定2**:

$$
|\mathcal R_i(t,x)| \le a\,\|x\|_r^{\,r_i+\varrho}
    \qquad(\forall\,t+\tau\ge0,\ \forall x\in\mathbb R^n,\ i=1,\dots,n)
$$

<div align="right">式 (5-2)</div>

を満たすものとする($a>0$, $\varrho+r_{\min}>0$; これは各時刻 $t$ を固定するごとに $\mathcal R(t,\cdot)$ が次数 $\varrho$ の $r$-同次ベクトル場であれば自動的に成り立つ評価である)。

**定理2**: 仮定1・仮定2のもとで, (i) $\nu>0$ かつ $\varrho>\nu$ ならば式(5-1)のゼロ解はAS。(ii) $\nu<0$ かつ $\varrho<\nu$ ならば式(5-1)の解はUUB。

## 1.2 貯留関数法の基礎式


#### 星・山岡モデルによる分布遅延型拡張

星・山岡(1982)は非定常項 $K_2\,\mathrm{d}(Q^{P_2})/\mathrm{d}t$ を追加した

$$
S(t) = K_1\,Q(t)^{P_1} + K_2\,\frac{\mathrm{d}}{\mathrm{d}t}\bigl[Q(t)^{P_2}\bigr]
$$

<div align="right">式 (6)</div>

という式を示した(「＜参考資料３＞貯留関数法とその適用法」([kohyo-21-k133-1-3](https://www.scj.go.jp/ja/info/kohyo/kohyo-21-k133.html)), pp.33–35 による)。ただし，$K_1,K_2,P_1,P_2$ は流域(または河道)パラメータ, $S$：貯留量, $Q$：流出量(または直接流出量) である。

式(6)を, 角屋・永井(1978)の水収支式を現在時刻 $t$ 基準・移動窓 $[t-\tau,t]$ の形に書き換えた式として，

$$
S(t) = \int_{t-\tau}^{t} I(s)\,\mathrm{d}s - \int_{t-\tau}^{t} Q(s)\,\mathrm{d}s
$$

<div align="right">式 (7)</div>

を得る。式(6)と式(7)より

$$
K_1\,Q(t)^{P_1} + K_2\,\frac{\mathrm{d}}{\mathrm{d}t}\bigl[Q(t)^{P_2}\bigr]
    = \int_{t-\tau}^{t} I(s)\,\mathrm{d}s - \int_{t-\tau}^{t} Q(s)\,\mathrm{d}s
$$

<div align="right">式 (8)</div>

が得られる。状態変数を次のように設定すれば，

$$
X(t) := Q(t)^{P_2}, \qquad Q(t) = X(t)^{1/P_2}
$$

<div align="right">式 (9)</div>

$\mathrm{d}X/\mathrm{d}t = \mathrm{d}(Q^{P_2})/\mathrm{d}t$ であるから, 式(8)は $X$ についての連続時間状態方程式である

$$
\dot X(t) = -\frac{K_1}{K_2}\,X(t)^{P_1/P_2}
    + \frac1{K_2}\int_{t-\tau}^{t} I(s)\,\mathrm{d}s
    - \frac1{K_2}\int_{t-\tau}^{t} X(s)^{1/P_2}\,\mathrm{d}s
$$

<div align="right">式 (10)</div>

が得られる。

$$
F(X) := -\frac{K_1}{K_2}\,X^{P_1/P_2}, \qquad
  R(t) := \frac1{K_2}\int_{t-\tau}^{t}\Bigl[I(s)-X(s)^{1/P_2}\Bigr]\,\mathrm{d}s
$$

<div align="right">式 (11)</div>

とおくと $\dot X(t)=F(X(t))+R(t)$。これは§1.1の定理2(式(5-1))において $G\equiv0$ とし, 摂動 $\mathcal R$ を $R(t)$ の被積分関数とみなした場合に対応する。同定理は, 摂動 $\mathcal R$ が仮定2(式(5-2))を満たすとき, $\nu>0,\varrho>\nu$ でゼロ解 AS, $\nu<0,\varrho<\nu$ で解UUBとなる。

#### 仮定1の検証(本モデルへの適用)

本モデルには $\int_{t-\tau}^t G(x(s))\,\mathrm{d}s$ に相当する項が存在しない ($R(t)$ の全体を摂動として扱うため $G\equiv0$)。したがって §1.1の仮定1が要求する補助系(3)は $G\equiv0$ の場合として $\dot X=F(X)$ 自体に一致し, 式(11) より

$$
\dot X(t) = -\frac{K_1}{K_2}\,X(t)^{P_1/P_2}
$$

<div align="right">式 (12)</div>

となる。係数 $-K_1/K_2<0$ は $K_1,K_2>0$(§1.2 パラメータの制約)より任意の $K_1,K_2$ に対し成り立つから, ゼロ解は大域的AS (基本状態方程式 $\dot Q=-Q^{1/p}$ に対する対応する補助系の議論と同型)。よって仮定1は追加条件なしに自動的に成立する。

#### パラメータの制約

「＜参考資料３＞貯留関数法とその適用法」([kohyo-21-k133-1-3](https://www.scj.go.jp/ja/info/kohyo/kohyo-21-k133.html), pp.33–35)によれば, 式(6)のパラメータには次の制約がある。

- $K_1,K_2>0$, $P_1,P_2>0$

<div align="right">式 (13)</div>

#### 仮定2の検証(本モデルへの適用)

- **主項の同次次数の符号**: 式(11)の主項の次数は $\nu=P_1/P_2-1$ である。 $\dot X=F(X)$ が定理1(i)の意味でAS favorable($\nu>0$)となる条件は

$$
P_1 > P_2
$$

<div align="right">式 (14)</div>

それに加え，$R(t)$ の被積分関数のうち $X(s)^{1/P_2}$ の次数は $\varrho=1/P_2-1$ である。降雨 $I(s)$ は, 洪水イベント全体での流出係数(総流出量/総降雨量)がおおむね1に近いという経験則により, $I(s)$ は $Q(s)=X(s)^{1/P_2}$ と同程度の大きさをもつとみなせる。 定理2のAS条件は $\varrho>\nu$ ($\nu=P_1/P_2-1$)という**狭義**の不等号を貯留係数に当てはめると
$$
\frac1{P_2}-1 > \frac{P_1}{P_2}-1 \iff P_1<1
$$

<div align="right">式 (15)</div>

である。上記条件は貯留の制約が($P_1<1$)であるため, 式(14)($P_1>P_2$, $\nu>0$)とあわせて, 定理2によりゼロ解ASが成立する。

- **文献の代表値との整合**: 星・村上(1987)の代表値 $P_1=0.6,P_2=0.4648$(式(13))は $0.6>0.4648$ より条件(14)を満たす。

## 1.3 MRAC型適応制御の式展開

問題設定は

$$
\dot X(t) = -\theta_3\,X(t)^{\alpha_1} + R(t)
$$

<div align="right">式 (16)</div>

$$
R(t) := \theta_1\Phi_I(t) - \theta_2\Phi_Q(t)
$$

<div align="right">式 (17)</div>

$$
\Phi_I(t) := \int_{t-\tau}^{t} I(s)\;=\;\int_{t-\tau}^{t} X(s)^{\alpha_2}\,\mathrm{d}s
$$

<div align="right">式 (18)</div>

$$
\Phi_Q(t) := \int_{t-\tau}^{t} X(s)^{\alpha_2}\,\mathrm{d}s
$$

<div align="right">式 (19)</div>

とする。式(16)は§1.2の $\dot X=F(X)+R(t)$(式(10、11))と同形であり, $F(X):=-\theta_3X^{\alpha_1}$ とおけば,  $K_1/K_2\to\theta_3$ の対応でそのまま流用できる。

$\theta_2=1/K_2$ は $\theta_2\ll1$ となるため, 定理2の摂動として自然に減衰する。したがって $\theta_1,\theta_2$ を個別に較正・推定する必要はなく, $\theta_3$ のみをMRACで追跡すれば十分である。

#### 誤差方程式

状態誤差とパラメータ誤差を

$$
e(t):=X(t)-\hat X(t)
$$

<div align="right">式 (23)</div>

$$
\tilde\theta_3:=\hat\theta_3-\theta_3
$$

<div align="right">式 (24)</div>

($e(t)$\[m$^3$/s\] は状態誤差, $\tilde\theta_3$はパラメータ誤差)と定義すると, 式(16)から(22)を引くと, 恒等的に相殺し

$$
\dot e = -\bigl(\theta_3 X^{\alpha_1}-\hat\theta_3\hat X^{\alpha_1}\bigr)
$$

<div align="right">式 (25)</div>

が成り立つ。$\theta_3=\hat\theta_3-\tilde\theta_3$ を代入して整理すると

$$
\dot e = -\hat\theta_3\bigl(X^{\alpha_1}-\hat X^{\alpha_1}\bigr) + \tilde\theta_3 X^{\alpha_1}
$$

<div align="right">式 (26)</div>

を得る。

#### Lyapunov関数と適応則

適応ゲイン $\gamma_3>0$(無次元)を導入し, 候補Lyapunov関数を

$$
V = \tfrac12 e^2 + \tfrac{1}{2\gamma_3}\tilde\theta_3^2
$$

<div align="right">式 (28)</div>

とおくと, 式(26)より

$$
\dot V = -\hat\theta_3\,e\bigl(X^{\alpha_1}-\hat X^{\alpha_1}\bigr)
    + \tilde\theta_3\Bigl[\,e\,X^{\alpha_1}+\tfrac{1}{\gamma_3}\dot{\hat\theta}_3\Bigr].
$$

<div align="right">式 (29)</div>

$[\,\cdot\,]$ 内を恒等的にゼロにする適応則を選ぶことで交差項を相殺する:

$$
\dot{\hat\theta}_3(t) = -\gamma_3\,e(t)\,X(t)^{\alpha_1}
$$

<div align="right">式 (30)</div>


#### 遅延積分項とLyapunov–Krasovskii汎関数について

$R(t)$(したがって $\Phi_I(t),\Phi_Q(t)$)は分布遅延の窓積分を含むが, 式(25)で見た通り真の系・推定器の双方に同一の値で現れ誤差方程式 $\dot e$ から完全に消去されるため, 本節の $V=\tfrac12e^2+\tfrac1{2\gamma_3}\tilde\theta_3^2$ による解析は遅延項の影響を一切受けない。

# 2. sinc カーネル遅延埋め込みレプリケータ動力学によるピーク遅れ時間推定
## 2.2 レプリケータ関数

 「候補遅延 $k\Delta\tau$ だけ過去に遡った降雨 $I(t-k\Delta\tau)$ が, 現在の流出量 $\hat Q(t)$ にどれだけsinc類似度で一致するか」を各タップ $k$ の適応度とする, 直接照合型のレプリケータ選択に組み替える。
両信号をそれぞれ自身のピーク値で正規化した無次元量とする。

$$
I_{\mathrm{norm}}(t) := \frac{I(t)}{\max_t I(t)}
$$

<div align="right">式 (41)</div>

$$
Q_{\mathrm{norm}}(t) := \frac{\hat Q(t)}{\max_t \hat Q(t)}
$$

<div align="right">式 (42)</div>

候補遅延 $M$ 個($\Delta\tau$刻み, $k=0,\dots,M-1$)による降雨側の遅延埋め込みを

$$
 x(t)[k] := I_{\mathrm{norm}}(t-k\Delta\tau)
$$

<div align="right">式 (43)</div>

とし, タップ $k$ の適応度を, sinc カーネルを用いた 

$$
f_k(t) := \mathrm{sinc}\bigl(x(t)[k] - Q_{\mathrm{norm}}(t)\bigr)
$$

<div align="right">式 (44)</div>

で定義する。

$$
\dot\omega_k(t) = \omega_k(t)\Bigl(f_k(t) - \omega(t)^{\!\top}f(t)\Bigr)
$$

<div align="right">式 (45)</div>

$$
\omega_k(0) = \frac1M \qquad(k=0,\dots,M-1)
$$

<div align="right">式 (46)</div>

とする。式(45)を時刻 $0$ から観測終了時刻 $T$ まで数値積分し, 終端の重み $\omega(T)$ が最大となる値を

$$
k^\ast := \mathop{\mathrm{argmax}}_k \omega_k(T)
$$

<div align="right">式 (47)</div>

$$
t_{\mathrm{lag}} := k^\ast\,\Delta\tau
$$

<div align="right">式 (48)</div>

を推定ピーク遅れ時間とする。$k^\ast$ が候補遅延の最大タップ($M-1$)に一致した場合は, 真の遅れがカバー範囲を超えているため, タップ数 $M$ を増やして再実行すべき旨を警告する。

## 2.3 次回移動窓長 $\tau$ の目安式

推定されたピーク遅れ時間 $t_{\mathrm{lag}}$ と, 既知の貯留定数 $K_{\min}$(`CatchmentModel.K_min`)から, 次回のシミュレーションで設定し直す際の目安として

$$
\tau_{\mathrm{suggest}} := 2\bigl(t_{\mathrm{lag}} - K_{\min}\bigr)
$$

<div align="right">式 (49)</div>

を計算・出力する。

#### 目安式の導出

以下は式(49)の近似的な導出である。本モデルの非線形性($Q^{1/p}$ 項)を無視した線形近似のもとでの見積であり, 厳密解ではないことに注意する。§6.1の基本状態方程式(式(52))を, 概念的に次の2段階の直列系とみなす。

**第1段階(移動窓平滑化)**: $w(t) := \frac1{K_2}\int_{t-\tau}^{t} I(s)^{1/p}\,\mathrm{d}s$ は, 幅 $\tau$ の矩形移動平均フィルタ(インパルス応答 $h_w(u)=1/\tau\ (0\le u\le\tau)$, $u$:遅れ)とみなせる。このインパルス応答の重心(平均遅れ時間)は

$$
\bar t_w = \int_0^\tau u\cdot\frac{1}{\tau}\,\mathrm{d}u = \frac{\tau}{2}
$$

<div align="right">式 (49-1)</div>

**第2段階(貯留・ルーティング)**: 式(52)の窓積分減衰項 $-\frac1K\int Q(s)^{1/p}\mathrm{d}s$ と瞬時減衰項 $-Q^{1/p}$ を, 貯留定数 $K_{\min}$ をもつ実効的な1個の線形貯留系(インパルス応答 $h_K(u)=\frac1{K_{\min}} e^{-u/K_{\min}}\ (u\ge0)$)にまとめて近似する。このインパルス応答の平均滞留時間は

$$
\bar t_K = \int_0^\infty u\cdot\frac1{K_{\min}} e^{-u/K_{\min}}\,\mathrm{d}u = K_{\min}
$$

<div align="right">式 (49-2)</div>

直列につながる線形系では, 全体の遅れ時間は各段階の遅れ時間の和になる。したがって

$$
t_{\mathrm{lag}} \;\approx\; \bar t_w + \bar t_K \;=\; \frac{\tau}{2} + K_{\min}
$$

<div align="right">式 (49-3)</div>

を $\tau$ について解くと, 式(49)が得られる。

$\tau_{\mathrm{suggest}}$ は移動窓の*長さ*であり正でなければならないため, 式(49)は

$$
t_{\mathrm{lag}} - K_{\min} > 0 \qquad\text{すなわち}\qquad t_{\mathrm{lag}} > K_{\min}
$$

<div align="right">式 (50)</div>

とする。$t_{\mathrm{lag}}\le K_{\min}$(推定遅れ時間が貯留定数以下)の場合は $\tau_{\mathrm{suggest}}\le0$ となり式の適用範囲外であるため, この目安式は採用せず, 遅延推定自体(タップ数 $M$ や $\Delta\tau$ の設定)を見直すべきである。採用するかどうかは, 出力されるグラフ・数値(降雨・流出量波形, レプリケータ重みの収束, sinc適応度の推移)を目視で確認した上で判断する必要がある。

# 3. 全体構成

パイプライン全体は次の順序で構成される。

1.  **降雨波形の生成**: 三角形型のハイエトグラフ $I(t)$ を合成する。

2.  **流出量のオンライン推定**: 分布遅延をもつ貯留関数型の状態方程式 $\dot X(t)$(式(16), §1.3)に対し, 未知パラメータ $\theta_3$ のみをMRAC型適応則で推定する並列モデル $\hat X(t)$ (式(21))を数値積分する。

3.  **ポンプ場への流入**: $\hat Q(t)$ を1秒刻みに補間し, 可変速ポンプ+ ウェットウェルからなるポンプ場モデルへの流入量として与え, 定速(CS) ヒステリシス制御下での吐出量 $Q_R(t)$・消費電力を計算する。

4.  **ピーク遅れ時間の推定**: $I(t)$ と $\hat Q(t)$ (5分刻み) から, sinc カーネル遅延埋め込み + レプリケータ動力学により両者のピーク遅れ時間 $t_{\mathrm{lag}}$ を推定する。

現在の実装は, 手順3のポンプ吐出量 $Q_R(t)$ は手順2の $\hat Q(t)$ の計算には*フィードバックされない* $\theta_1\Phi_I$ 項に反映させる拡張は実装上は存在するが, 現在のパイプラインでは無効化されている。

すなわち全イベント・全時刻を通じて $I_{\mathrm{flow}}(t)\in[64,\ 4051]$m$^3$/min (非零区間; m$^3$/s換算でおよそ1–68m$^3$/s)のスケールに収まる。この $I_{\mathrm{flow}}(t)$ の桁(alpha2乗されて回帰信号 $\Phi_I$ に入る際にさらに増幅される)に応じて, $\theta_3$ の較正値も大きく変わる(本稿執筆時点の実測値は `PipelineConfig.adaptive_theta3_true`$\approx1.0\times10^7$)。

# 4. 記法

<div class="center">

| 記号                       | 単位                  | 意味                                                                                                                                                                                                                                          |
|:---------------------------|:----------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| $I(t)$                     | mm                    | 降雨強度(観測値, $\Delta t_{\min}$刻みの積算降雨量)                                                                                                                                                                                           |
| $Q(t)$                     | m$^3$/s               | 流出量(真値)                                                                                                                                                                                                                                  |
| $\hat Q(t)$                | m$^3$/s               | 流出量の推定値, $:=\hat X(t)^{1/P_2}$(式(22))                                                                                  |
| $X(t)$                     | (m$^3$/s)$^{P_2}$     | §1.3の状態, $:=Q(t)^{P_2}$(式(9)と同じ変換; 物理量そのものではない) |
| $\hat X(t)$                | (m$^3$/s)$^{P_2}$     | $X(t)$ の並列モデル推定値(式(21)で積分される状態そのもの)                                                                                                       |
| $p$                        | –                     | 積分項の指数を決める形状パラメータ(§6.1), $\alpha:=1/p$                                                                                                     |
| $\alpha_1$                 | –                     | 瞬時減衰項の指数, $:=P_1/P_2$(既知・固定)                                                                                                                                                                                                     |
| $\alpha_2$                 | –                     | 窓積分項の指数, $:=1/P_2$(既知・固定)                                                                                                                                                                                                         |
| $A$                        | km$^2$                | 流域面積(既知)                                                                                                                                                                                                                                |
| $\theta_1$                 | 1/min                 | 流入寄与係数 $A^{\alpha_2}/K_2$(既知の物理定数, 推定しない)                                                                                                                                                                                   |
| $\theta_2$                 | 1/min                 | 減衰寄与係数 $1/K_2$(既知の物理定数, 推定しない; $\theta_2\ll1$)                                                                                                                                                                              |
| $\theta_3,\hat\theta_3(t)$ | –                     | 瞬時減衰項 $-\theta_3 X^{\alpha_1}$ の係数 $K_3$(真値, 推定値; 唯一の未知パラメータ)                                                                                                                                                          |
| $\tau$                     | min                   | 分布遅延の積分窓長(移動窓長)                                                                                                                                                                                                                  |
| $R(t)$                     | (m$^3$/s)$^{P_2}$/min | 既知の外生入力 $\theta_1\Phi_I(t)-\theta_2\Phi_Q(t)$(式(17))                                                                               |
| $\Phi_I(t),\Phi_Q(t)$      | –                     | 回帰信号(降雨側・流出量側の窓積分, $\Phi_Q$ の被積分関数は式 (19)より $X(s)^{\alpha_2}=Q(s)$)                                                    |
| $e(t)$                     | (m$^3$/s)$^{P_2}$     | 状態誤差 $X(t)-\hat X(t)$(物理量Q(t)の誤差ではない, 式(23))                                                                                            |
| $\gamma_3$                 | –                     | 適応ゲイン                                                                                                                                                                                                                                    |
| $m(t)^2$                   | –                     | 正規化勾配法の共通分母(未使用)                                                                                                                                                                                                                |
| $A$                        | km$^2$                | 対象流域面積                                                                                                                                                                                                                                  |
| $c_A$                      | m$^3$/(mm$\cdot$min)  | 降雨深\[mm\]$\to$流入水量\[m$^3$/s\]換算係数                                                                                                                                                                                                  |
| $K_{\min}$                 | min                   | 貯留関数法の貯留定数(`CatchmentModel.K_min`)                                                                                                                                                                                                  |
| $N$                        | rpm                   | ポンプ回転数                                                                                                                                                                                                                                  |
| $Q_{\mathrm{op}}$          | m$^3$/s               | ポンプ動作点流量(揚程曲線とプラント曲線の交点)                                                                                                                                                                                                |
| $H_N(Q_{\mathrm{op}})$     | m                     | 回転数 $N$ における揚程曲線                                                                                                                                                                                                                   |
| $\eta_N(Q_{\mathrm{op}})$  | –                     | 回転数 $N$ における効率                                                                                                                                                                                                                       |
| $P$                        | W                     | ポンプ消費電力                                                                                                                                                                                                                                |
| $H_w(t)$                   | m                     | ウェットウェル水位                                                                                                                                                                                                                            |
| $S$                        | m$^2$                 | ウェットウェル断面積                                                                                                                                                                                                                          |
| $Q_R(t)$                   | m$^3$/s               | ポンプ吐出量(リアルタイム結合における流出側)                                                                                                                                                                                                  |
| $\omega_k(t)$              | –                     | レプリケータ重み(候補遅延 $k$ タップ目のシンプレックス成分)                                                                                                                                                                                   |
| $\Delta\tau$               | min                   | 候補遅延タップの間隔                                                                                                                                                                                                                          |
| $M$                        | –                     | 候補遅延タップ数                                                                                                                                                                                                                              |
| $t_{\mathrm{lag}}$         | min                   | 推定ピーク遅れ時間                                                                                                                                                                                                                            |

</div>

$\tau$(貯留状態方程式の移動窓長)と $\Delta\tau$(遅延タップ間隔)は実装上どちらも `tau_min` という変数名で登場するが別物であるため, 本稿では記号を分けて区別する。またポンプ節の流量 $Q_{\mathrm{op}}$ は「ポンプの動作点における吐出流量」であり, §6 の流出量 $Q(t)$ とは物理的に別の量である(実装コードではいずれも変数名 `Q` が使われている)。

# 5. 降雨波形モデル

出典: 横浜市環境創造局,「横浜市下水道計画指針-2010年版-」[第6章 雨水管理計画](https://www.city.yokohama.lg.jp/kurashi/machizukuri-kankyo/kasen-gesuido/gesuido/shishin_2010.files/shishin2010-06.pdf), §6.3.3 確率年別降雨強度式・§6.4.2 ハイエトグラフ(目次ページ: [横浜市下水道計画指針-2010年版-](https://www.city.yokohama.lg.jp/kurashi/machizukuri-kankyo/kasen-gesuido/gesuido/shishin_2010.html))。本節の三角形型ハイエトグラフは, 同指針§6.4.2が採用する中央集中型ハイエトグラフ(図6.4.2.2)を模したものである。§1.3の $I_{\mathrm{flow}}(t)$ のスケール検証(表)で用いた5/10/20/30年確率のIDF式(降雨強度式)も同指針§6.3.3(表6.3.3.1)による。

離散時刻 $t_k = k\,\Delta t_{\min}$ ($k=0,\dots,n-1$, $\Delta t_{\min}$は観測周期, 既定5分) 上で, ピーク時刻 $t_{\mathrm{peak}}$, 立ち上がり時間 $T_{\mathrm r}$, 立ち下がり時間 $T_{\mathrm f}$, ピーク強度 $I_{\mathrm{peak}}$ をもつ三角形型ハイエトグラフを

$$
I(t_k) =
  \begin{cases}
    \displaystyle I_{\mathrm{peak}}\,\frac{t_k-(t_{\mathrm{peak}}-T_{\mathrm r})}{T_{\mathrm r}}
      & t_{\mathrm{peak}}-T_{\mathrm r}\le t_k \le t_{\mathrm{peak}} \\[4pt]
    \displaystyle I_{\mathrm{peak}}\,\frac{(t_{\mathrm{peak}}+T_{\mathrm f})-t_k}{T_{\mathrm f}}
      & t_{\mathrm{peak}} < t_k \le t_{\mathrm{peak}}+T_{\mathrm f} \\[4pt]
    0 & \text{それ以外}
  \end{cases}
$$

<div align="right">式 (51)</div>

により生成し, 負値は0にクリップする。乱数シードを指定した場合は, $I(t_k)>0$ の各点に独立な加法的ガウス雑音 $\mathcal N(0,(0.05\,I_{\mathrm{peak}})^2)$ を加え, 再度0にクリップする。既定設定では小雨イベント($I_{\mathrm{peak}}=5$mm)・大雨イベント($I_{\mathrm{peak}}=20$mm)の2ケースを, 同一の $t_{\mathrm{peak}}=90$分, $T_{\mathrm r}=50$分, $T_{\mathrm f}=70$分, 継続時間240分で生成する。

# 6. 流出量モデル: 分布遅延をもつ貯留関数系のオンライン推定

## 6.1 基本状態方程式

対応: `appendix_q_model_simulation.simulate_appendix_q_model`

流出量 $Q(t)$ の基礎モデルは,状態方程式

$$
\dot Q(t) = -Q(t)^{1/p} + \frac{1}{K_2}\int_{t-\tau}^{t} I(s)^{1/p}\,\mathrm{d}s
                       - \frac{1}{K}\int_{t-\tau}^{t} Q(s)^{1/p}\,\mathrm{d}s
$$

<div align="right">式 (52)</div>

で与えられる。$p>0$は形状パラメータで, $p<1$(同次次数 $\nu=1/p-1>0$) のとき零解は漸近安定であることが, 同次システムの安定性に関する定理1(i) (§1.1)により示されている。 $K$\[min\] は貯留定数, $\tau$\[min\] は分布遅延の幅である。$K_2$\[min\] は

$$
K_2 := \frac{K}{A\cdot 1000/60}
$$

<div align="right">式 (53)</div>

$$
\frac{1}{K_2} = A\cdot\frac{1000}{60}\cdot\frac{1}{K}
$$

<div align="right">式 (54)</div>

で定義され, $A$\[km$^2$\]は対象流域面積である。係数 $1000/60$ は降雨深\[mm\]$\times$ 流域面積\[km$^2$\] $\to$ 流入水量\[m$^3$/s\]の単位換算(降雨深1mm$\times$面積1km$^2$ $=$ 水量1000m$^3$, 観測周期の分$\to$秒換算60)である。以下, この換算係数を

$$
c_A := A\cdot\frac{1000}{60}
$$

<div align="right">式 (55)</div>

と書く。既定では  流域面積 $A=0.45+0.55+0.30+0.20=1.5$km$^2$ を用いる。

式(52)は硬い非線形の瞬時減衰項 $-Q^{1/p}$ を含むため, 前進Euler(陽解法)では $\Delta t_{\min}=5$分程度の粗い刻みで振動的な数値不安定 (棘状のオーバーシュート)を起こす。そこで瞬時減衰項のみを後退Euler(陰解法)で扱い, 分布遅延積分項は台形則により陽的に評価する:

$$
Q[t] + \Delta t_{\min}\,Q[t]^{1/p} = Q[t-1] + \Delta t_{\min}\left(
    \frac{1}{K_2}\Phi_I[t] - \frac{1}{K}\Phi_Q[t]\right)
$$

<div align="right">式 (56)</div>

を $Q[t]\ge 0$ について求根(`scipy.optimize.brentq`)する。ここで $\Phi_I[t],\Phi_Q[t]$ は時刻 $t-1$ までの履歴を用いた窓積分の台形則近似 (§6.2参照)である。

式(52)では瞬時減衰項の係数が暗に $1$ に固定されているが, これを未知パラメータ $\theta_3:=K_3$ として一般化し, 降雨側・流出側の係数 $\theta_1,\theta_2$ は既知の物理定数として扱った上でMRAC型オンライン推定則(真の系・並列モデル, 誤差方程式, Lyapunov関数, 適応則)により $\theta_3$ のみを推定する式展開は§1.3(第1章)にまとめてある。以下の $\Phi_I,\Phi_Q$ 等の記法はそこでの定義(式(17), 指数 $\alpha_2:=1/P_2$)をそのまま用いる。

## 6.2 回帰信号の離散化

窓長 $\tau$ を $\Delta t_{\min}$ 刻みで離散化した点数 $w:=\mathrm{round}(\tau/\Delta t_{\min})$ を用い, 時刻$t$ における回帰信号は, 開始点 $s_0:=\max(0,t-w)$ として台形則(`numpy.trapz`)により

$$
\Phi_I[t] \approx \int_{s_0}^{t-1} I[u]^{\alpha_2}\,\mathrm{d}u
$$

<div align="right">式 (57)</div>

$$
\Phi_Q[t] \approx \int_{s_0}^{t-1} Q[u]^{\alpha_2}\,\mathrm{d}u
$$

<div align="right">式 (58)</div>

と近似する。現在時刻 $t$ の値はまだ確定していないため, 瞬時項の後退Euler更新 (§6.1)と整合させて, 履歴は $t-1$ までで打ち切る。

## 6.3 数値積分アルゴリズムの要約

真の系(式(16))・推定器(式(21))はいずれも, 瞬時減衰項のみ後退Euler+分布遅延積分項は台形則による陽的評価, という半陰的Euler法で数値積分する。適応則(式(30))は前進Eulerで積分する。3つのモード

- `"true"`: 真値 $\theta_1,\theta_2,\theta_3$($\theta_1,\theta_2$ は既知定数, $\theta_3$ も真値)で積分した参照軌道そのもの(推定器なし)。

- `"fixed"`: $\hat\theta_3$ を誤った初期値(真値の比 `*_hat0_ratio` 倍)に固定し, 適応しない($\gamma_3=0$)。 $\theta_1,\theta_2$ は真の系と同じ既知定数を用いる。

- `"adaptive"`: `"fixed"`と同じ初期値から, 式 (30)で $\gamma_3$ を用いてオンライン推定する。

を同一の設定から生成できる。


# 7. ポンプ場(ウェットウェル)モデル

 ※※(Fecarotta et al., 2018, Resources 7, 73)より丸々引用。台数を増やした拡張 (§7.6)以外は, 必要最小限の条件・パラメータ式のみ引用する。

## 7.1 ポンプ性能曲線

揚程曲線($N$\[rpm\]:回転数, $Q_{\mathrm{op}}$\[m$^3$/s\]:吐出流量): 

$$
H_N(Q_{\mathrm{op}}) = c_{0h}N^2 + c_{2h}Q_{\mathrm{op}}^2
$$

<div align="right">式 (60)</div>

 締切揚程比 $r_{\mathrm{so}}$(無次元)とBEP点(最高効率点; $Q_{\mathrm{op}}=
Q_{\mathrm{bep}}(N_{\max})$, $H=H_{\mathrm{bep}}(N_{\max})$)の2条件により: 

$$
c_{0h} = \frac{r_{\mathrm{so}}H_{\mathrm{bep}}(N_{\max})}{N_{\max}^2}
$$

<div align="right">式 (61)</div>

$$
c_{2h} = \frac{(1-r_{\mathrm{so}})H_{\mathrm{bep}}(N_{\max})}{Q_{\mathrm{bep}}(N_{\max})^2}
$$

<div align="right">式 (62)</div>

 BEP効率($r_\eta$\[無次元\]:$N_{\min}$での低下率): 

$$
\eta_{\mathrm{bep}}(N) = \eta_{\mathrm{bep}}(N_{\max}) - \frac{r_\eta\,\eta_{\mathrm{bep}}(N_{\max})}{N_{\max}-N_{\min}}(N_{\max}-N)
$$

<div align="right">式 (63)</div>

 相対効率($x:=Q_{\mathrm{op}}/N$, $x_{\mathrm{bep}}:=Q_{\mathrm{bep}}(N_{\max})/N_{\max}$, $c_{\mathrm{curv}}>0$:曲率, $e_{\mathrm{floor}}\in(0,1)$:下限, いずれも無次元): 

$$
e(x) = \mathrm{clip}\!\left(1 - c_{\mathrm{curv}}\left(\frac{x-x_{\mathrm{bep}}}{x_{\mathrm{bep}}}\right)^{2},\;
    e_{\mathrm{floor}},\,1\right)
$$

<div align="right">式 (64)</div>

 効率と消費電力($\gamma_w=9806$N/m$^3$:水の比重量, $P$\[W\]): 

$$
\eta_N(Q_{\mathrm{op}}) = e(Q_{\mathrm{op}}/N)\,\eta_{\mathrm{bep}}(N)
$$

<div align="right">式 (65)</div>

$$
P(Q_{\mathrm{op}},N) = \frac{\gamma_w\,H_N(Q_{\mathrm{op}})\,Q_{\mathrm{op}}}{\eta_N(Q_{\mathrm{op}})}
$$

<div align="right">式 (66)</div>

## 7.2 プラント諸量

$\beta:=H_0/H_{\mathrm{bep}}(N_{\max})$(無次元)により静水頭 $H_0$\[m\]と配管損失係数 $K$\[m$\cdot$s$^2$/m$^6$\]を定める: 

$$
H_0 = \beta\,H_{\mathrm{bep}}(N_{\max})
$$

<div align="right">式 (67)</div>

$$
K = \frac{H_{\mathrm{bep}}(N_{\max})-H_0}{Q_{\mathrm{bep}}(N_{\max})^2}
$$

<div align="right">式 (68)</div>

 $S_{h,\max}$\[回/h\]:1時間あたり最大起動回数, $S$\[m$^2$\]:ウェットウェル断面積($H_{w,\min}=0$)として: 

$$
H_{w,\max} = \frac{900\,Q_{\mathrm{bep}}(N_{\max})}{S_{h,\max}\,S}
$$

<div align="right">式 (69)</div>

 目標水位 $H_{w,\max}^{\mathrm{target}}$\[m\]からの逆算(`build_pump_station`): 

$$
S = \frac{900\,Q_{\mathrm{bep}}(N_{\max})}{S_{h,\max}\,H_{w,\max}^{\mathrm{target}}}
$$

<div align="right">式 (70)</div>

 最大許容吐出量($r_{Qad}$:無次元比例定数): 

$$
Q_{\mathrm{ad,max}} = r_{Qad}\,Q_{\mathrm{bep}}(N_{\max})
$$

<div align="right">式 (71)</div>

## 7.3 動作点方程式

揚程曲線(式(60))とプラント曲線の交点(交点は高々1つ; 実行不可能なら $Q_{\mathrm{op}}$ は不定, 上限 $Q_{\mathrm{ad,max}}$ でクリップ): 

$$
H_0-H_w+K\,Q_{\mathrm{op}}^2 = c_{0h}N^2+c_{2h}Q_{\mathrm{op}}^2
$$

<div align="right">式 (72)</div>

$$
Q_{\mathrm{op}}^2 = \frac{H_0-H_w-c_{0h}N^2}{c_{2h}-K}
$$

<div align="right">式 (73)</div>

## 7.4 定速(CS)ヒステリシス制御と水位更新

$$
\text{ON}\to\text{OFF} \quad\text{if}\quad H_w\le H_{w,\min}
$$

<div align="right">式 (74)</div>

$$
\text{OFF}\to\text{ON} \quad\text{if}\quad H_w\ge H_{w,\max}
$$

<div align="right">式 (75)</div>

 ON中は $N=N_{\max}$ で式(73), OFF中は $Q_{\mathrm{op}}=0,N=0$。水位はEuler更新(刻み $\Delta t_s$, 既定1秒): 

$$
H_w[i] = \mathrm{clip}\!\left(H_w[i-1] + \frac{Q_{\mathrm{in}}[i]-Q_{\mathrm{op}}[i]}{S}\,\Delta t_s,\;
    H_{w,\min},\,H_{w,\max}\right)
$$

<div align="right">式 (76)</div>

## 7.5 エネルギー指標

消費電力量とその理論下限: 

$$
E = \sum_i P_i\,\frac{\Delta t_s}{3.6\times10^6}\;\text{[kWh]}
$$

<div align="right">式 (77)</div>

$$
E_{\mathrm{ref}} = \int \gamma_w\,Q_{\mathrm{in}}(t)\bigl(H_0+K\,Q_{\mathrm{in}}(t)^2\bigr)\,
    \frac{\mathrm{d}t}{3.6\times10^6}\;\text{[kWh]}
$$

<div align="right">式 (78)</div>

## 7.6 複数台ポンプ(段階起動)への拡張


同一機種のポンプが $n_p$ 台あり, 水位上昇に応じて1台ずつ段階的に追加起動し, 稼働中の $m$ 台($0\le m\le n_p$)が共通のウェットウェル・共通の吐出管に並列運転する構成に一般化する。全台が同一機種・同一回転数 $N$ のため, 対称性より合計吐出量 $Q_{\mathrm{total}}$ を等分に分担する。

$$
Q_{\mathrm{each}} = \frac{Q_{\mathrm{total}}}{m}
$$

<div align="right">式 (79)</div>

 各台の揚程曲線に代入すると, $m$ 台分の揚程は 

$$
H = c_{0h}N^2 + c_{2h}\left(\frac{Q_{\mathrm{total}}}{m}\right)^{2}
$$

<div align="right">式 (80)</div>

 一方プラント側の揚程は台数によらず合計流量 $Q_{\mathrm{total}}$ で決まるため, 動作点方程式(72)は次のように一般化される: 

$$
H_0-H_w+K\,Q_{\mathrm{total}}^2 = c_{0h}N^2+c_{2h}\left(\frac{Q_{\mathrm{total}}}{m}\right)^{2}
$$

<div align="right">式 (81)</div>

$$
Q_{\mathrm{total}}^2 = \frac{H_0-H_w-c_{0h}N^2}{c_{2h}/m^2-K}
$$

<div align="right">式 (82)</div>

 ($m=1$ で式(72)(73)に厳密に一致し, `verify_multi_matches_single_pump` で数値的にも確認済み)。

台数の段階制御は, 水位範囲 $[H_{w,\min},H_{w,\max}]$ を $n_p$ 等分した $n_p+1$個の閾値 

$$
H_w^{(k)} = H_{w,\min} + \frac{k}{n_p}\bigl(H_{w,\max}-H_{w,\min}\bigr)
    \qquad(k=0,1,\dots,n_p)
$$

<div align="right">式 (83)</div>

 を用い, 前ステップの稼働台数 $m_{\mathrm{prev}}$・水位 $H_w^{\mathrm{prev}}$ から 

$$
m \leftarrow m_{\mathrm{prev}}; \quad
  H_w^{\mathrm{prev}}\ge H_w^{(m+1)} \text{ の間 } m\leftarrow m+1
$$

<div align="right">式 (84)</div>

$$
H_w^{\mathrm{prev}}\le H_w^{(m-1)} \text{ の間 } m\leftarrow m-1
$$

<div align="right">式 (85)</div>

 と決める($n_p=1$ で単一ポンプ版のON/OFF判定(式(74) (76))に一致)。水位更新・合計消費電力は式(76) と同型のまま $Q_{\mathrm{op}}\to Q_{\mathrm{total}}$ とし, 

$$
P_{\mathrm{total}} = m\cdot P\!\left(\frac{Q_{\mathrm{total}}}{m},\,N\right)
$$

<div align="right">式 (86)</div>

 ($m=0$ なら $Q_{\mathrm{total}}=P_{\mathrm{total}}=0$)で与える。

# 8. リアルタイム結合

 式 (22)により状態 $\hat X(t)$ , $\hat Q(t)=\hat X(t)^{1/P_2}$)で得た $\hat Q(t)$($\Delta t_{\min}=5$分刻み)を, (`resample_to_1s`)により1秒刻みの流入量 $Q_{\mathrm{in}}(t_{\mathrm{sec}})$ に変換し, §7のポンプ場モデルへ毎秒与える 。

$$
Q_R[i] = \frac{1}{|B_i|}\sum_{j\in B_i} Q_{\mathrm{op}}[j]
$$

<div align="right">式 (87)</div>

$$
B_i := \{\,j : t_{\mathrm{sec}}[j]\in[\,i\,\Delta t_{\min},\,(i+1)\Delta t_{\min})\,\}
$$

<div align="right">式 (88)</div>

($B_i$ は $i$ 番目の $\Delta t_{\min}$ 区間に属する1秒刻みインデックスの集合)。


# 9. まとめ

本リポジトリでは, 降雨波形の生成から, 分布遅延をもつ貯留関数型流出モデルの推定ができた。

# 10. えっせい
インフラ制御ついての願い
制御工学を専攻する学生の多くは、シミュレーションの対象として、コンベアなどの身近なものから、ロボット、ドローン、さらには月面探査機といった先端的な不安定物体を好んで設定する。しかし、コンベアを除けば、これらの多くは制御以外の分野に深刻な課題を抱えており、社会に広く普及しているとは言い難い。他方、インフラ分野に目を向けると、いまだに昭和時代の古い回転機器が現役で稼働している。維持管理も熟練の職人の経験に依存しており、自動化からはほど遠いのが現実である。インフラ現場空間（？）は先ほど言った先端機器制御の理想空間と全くちがう空間に存在しているものに見える。この「先端制御が描く理想空間」と「実際のインフラ現場」の乖離を埋めるため、現在インフラ界隈では制御システムの早急な刷新が求められている。インフラのメリットして、それを支える莫大な予算あるため、その一部を制御システムの開発費へと転換できれば、すぐさま自分の制御が現実のものとなる。だからこそ、、願わくは、制御学生はインフラを制御対象になりえるものとして頭の片隅に置いていてほしい。厚かましいようで申し訳ないのだが、一般化した制御系を開発したときには、シミュレーション対象にインフラを選んでいただければもう万々歳である。

# 参考文献

1.  山中理, 長岩明弘, 松原慎一郎, 仲田雅司郎, 山田富美夫, 「Hammerstein型非線形モデルを用いたシステム同定手法による下水道雨水流入量予測」, 電気学会論文誌D, Vol.120-D, No.4, pp.566-573, 2000. (§5の有効降雨算出・貯留関数法の基礎)

2.  O. Fecarotta, A. Carravetta, M.C. Morani, R. Padulano, “Optimal Pump Scheduling for Urban Drainage under Variable Flow Conditions,” *Resources*, 7, 73, 2018. doi:10.3390/resources7040073 (§7のポンプ性能曲線・ウェットウェルモデルの基礎)

3.  A. Aleksandrov, D. Efimov, E. Fridman, “Stability of homogeneous systems with distributed delay and time-varying perturbations,” *Automatica*, Vol. 153, 111058, 2023. (§1.1の定理1の出典)

4.   [横浜市下水道計画指針-2010年版-](https://www.city.yokohama.lg.jp/kurashi/machizukuri-kankyo/kasen-gesuido/gesuido/shishin_2010.html)

5.  星 清, 山岡 勲,「雨水流法と貯留関数法の相互関係」, 第26回水理講演会論文集, pp.273–278, 1982; 星 清, 村上泰啓,「小流域における総合貯留関数法の開発」, 第31回水理講演会論文集, pp.107–112, 1987 (いずれも「＜参考資料３＞貯留関数法とその適用法」(kohyo-21-k133-1-3), pp.33–35 に引用の式(20)–(24)を経由して参照; §1.2の星・山岡モデルによる分布遅延型拡張の基礎)

6.  Kim, W.H., Kim, Y.C., Ryu, J.W., “An Adaptive Storage Function Method for Rainfall-Runoff Forecasting,” *Transactions of the Korean Institute of Electrical Engineers*, Vol.47, No.2, pp.231–236, 1998.
    (§1のASFMの出典)

7.  Choi, S., Cho, T., Kim, W.H., Kim, Y.C., “An Adaptive Storage Function Method for Rainfall-Runoff Forecasting,” 計測自動制御学会論文集(Trans. of the Society of Instrument and Control Engineers), Vol.37, No.12, pp.1156–1161, 2001.
    (§1のASFMの出典、平昌川・忠州盆地の2流域での検証について)

---

## 免責事項(Disclaimer)

- 本リポジトリの内容(TeX/Markdown文書・数式・Pythonコード等)は，研究・学習目的で公開する解説資料であり，特定の流域・ポンプ場における実運用での有効性・安全性・正確性を保証するものではありません。
- 記載した数理モデル・安定性証明・パラメータの代表値は，引用文献に基づいて整理・実装したものです。個々の現場条件(流域特性，ポンプ諸元，制御盤仕様，安全基準等)への適合性は，利用者ご自身の責任で検証してください。
- 本リポジトリのコード・数式・記述内容を参考にプログラムを実装・改変・運用したことにより生じたいかなる損害(データ損失，誤動作，施設・設備への影響，人身・財産上の損害，金銭的損害等を含みますが，これらに限りません)についても，著作者は一切の責任を負いません。ご利用はすべて自己責任でお願いいたします。
- 本資料は法律・工学上の専門的助言を構成するものではありません。実際の下水道・排水ポンプ場等の設計・制御・運用に適用する場合は，必ず有資格の技術者・関係機関による検証・承認を経てください。
- 本リポジトリの著作物(文章・数式・図表・コード)の著作権は著者に帰属します。引用元として明記した文献(角屋・永井 1978；星・山岡 1982；星・村上 1987；Kim et al. 1998；Choi et al. 2001；Aleksandrov et al. 2023；Fecarotta et al. 2018 等)の著作権は各原著者・出版社に帰属し，本稿ではあくまで研究目的の範囲で要旨・数式を引用しています。第三者の著作物を含む再配布・転載を行う場合は，利用者の責任において原典の著作権者の許諾条件をご確認ください。

