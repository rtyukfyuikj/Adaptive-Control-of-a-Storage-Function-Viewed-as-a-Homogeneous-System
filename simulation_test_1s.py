
from __future__ import annotations

import importlib.util
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import matplotlib

# グラフをポップアップ表示するため, 利用可能なGUIバックエンドを探して指定する。
_GUI_BACKENDS = [
    ("qtagg", ("PyQt6", "PySide6", "PyQt5", "PySide2")),
    ("tkagg", ("tkinter",)),
]
for _backend, _required_any in _GUI_BACKENDS:
    if any(importlib.util.find_spec(_mod) is not None for _mod in _required_any):
        matplotlib.use(_backend)
        break
else:
    print("警告: PyQt/PySide/tkinterが見つからないため, グラフはポップアップ表示できず"
          "PNGファイルへの保存のみ行われます。")

import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

plt.rcParams["font.family"] = "MS Gothic"
plt.rcParams["axes.unicode_minus"] = False

UNIT_MM_KM2_TO_M3 = 1000.0
UNIT_MIN_TO_SEC = 60.0
GAMMA_WATER = 9806.0  # N/m^3, 水の比重量

# 出力PNGファイル名の接頭辞。スクリプトのファイル名から自動決定し, 他スクリプトとの
# 出力ファイル衝突を避ける。
OUTPUT_PREFIX = Path(__file__).resolve().stem

# np.trapzはNumPy 2.0でnp.trapezoidに改名されたため, 両バージョンに対応する。
trapezoid = getattr(np, "trapezoid", np.trapz)


# ----------------------------------------------------------------------
# 0. パラメータ集約
# 調整可能なパラメータはすべてこのdataclassにまとめてある。
# ----------------------------------------------------------------------
@dataclass
class PipelineConfig:
    # --- 対象流域 (rrl_storage_function_simulation.CatchmentModelにそのまま渡す) ---
    catchment_dt_min: float = 5.0
    catchment_time_area_km2: np.ndarray = field(
        default_factory=lambda: np.array([0.45, 0.55, 0.30, 0.20])
    )
    catchment_K_min: float = 15.0   # 修正RRL法+貯留関数法(Q_ref算出)の貯留定数K[分]

    # --- 降雨->流出量モデル (MRAC型オンライン適応推定) ---
    # theta1,theta2は既知の物理定数としてオンライン推定しない。推定するのはtheta3のみ。
    K_storage: float = 15.0   # theta1_true/theta2_true算出用の貯留定数K[分]
    p: float = 0.4              # P2: 窓積分項の指数 alpha2:=1/p (下のp1より小さく)
    p1: float = 0.60            # P1: 瞬時減衰項の指数 alpha1:=P1/P2 を決める形状パラメータ
    tau_min: float = 95.0     # 貯留状態方程式の移動窓長tau[分]


    adaptive_theta3_true: float = 1.0e7

    adaptive_theta3_hat0_ratio: float = 1.5
    adaptive_gamma3: float = 1.0e3
    adaptive_theta3_min: float = 1.0
    adaptive_theta3_max: float = 1.0e8

    # --- 降雨イベント (横浜市雨水管理計画 §6.3.3 確率年別降雨強度式)
    # I(t) = a/(t^n+b) [mm/hr], t=降雨継続時間[分]。5,10,20,30年確率の4式。
    rain_duration_min: float = 1440.0  # 降雨継続時間[分] (24時間, §6.4.2(1)の標準)
    rain_idf_params: dict = field(default_factory=lambda: {
        "5年確率": (880.0, 4.4, 0.65),
        "10年確率": (1452.0, 7.5, 0.70),
        "20年確率": (2199.0, 11.1, 0.75),
        "30年確率": (2731.0, 13.4, 0.77),
    })  # {ラベル: (a, b, n)} で I(t)=a/(t^n+b)

    # --- 降雨-流出量ピーク遅れ時間推定 (レプリケータ動力学+sincカーネル遅延埋め込み法) ---

    delay_tap_min: float = 5.0   # 候補遅延タップの間隔[分] (=catchment_dt_minに合わせる)
    n_taps: int = 19             # タップ数(候補遅延 0, delay_tap_min, ..., (n_taps-1)*delay_tap_min分)
    dt_eval: float = 1.0         # レプリケータODEの出力刻み[分]
    ivp_method: str = "RK23"
    ivp_max_step: float = 1.0
    ivp_rtol: float = 1e-6
    ivp_atol: float = 1e-8

    @property
    def delay_taps_min(self) -> np.ndarray:
        """候補遅延タップの遅延時間配列 [0, delay_tap_min, ..., (n_taps-1)*delay_tap_min]。"""
        return self.delay_tap_min * np.arange(self.n_taps)

    # --- ポンプ性能 (1台あたり, PumpParamsにそのまま渡す) ---

    pump_Nmax: float = 3000.0           # rpm (論文実測値)
    pump_Nmin: float = 1500.0           # rpm (論文実測値)
    pump_Q_bep_Nmax: float = 0.750       # m^3/s 
    pump_H_bep_Nmax: float = 12.0       # m (論文実測値)
    pump_H_shutoff_ratio: float = 10.0  # 限界流量
    pump_eta_bep_Nmax: float = 0.75
    pump_eta_bep_drop_at_Nmin: float = 0.15
    pump_rel_eff_curvature: float = 0.7
    pump_rel_eff_floor: float = 0.15

    # --- ウェットウェル/プラント設計 (PlantScenarioにそのまま渡す) ---
    plant_beta: float = 0.5
    plant_Sh_max: float = 10.0           # 1時間あたり最大起動回数 (論文実測値, 上記参照)
    plant_Hw_max_target: float = 100.0   # m (Hw_maxをこの値にするようSを逆算する)
    plant_Qad_max_ratio: float = 15.0

    # --- 複数台ポンプ構成(段階起動) ---
    n_pumps: int = 20

    # --- リアルタイム結合の時間刻み ---
    dt_s: float = 1.0   # [秒]


# ----------------------------------------------------------------------
# 1. 降雨生成 (横浜市雨水管理計画 §6.3.3降雨強度式 + §6.4.2中央集中型ハイエトグラフ)
#    + 修正RRL法+貯留関数法による参考流出量 (rrl_storage_function_simulation.py)
# ----------------------------------------------------------------------
def rainfall_intensity_idf(t_min: np.ndarray, a: float, b: float, n: float) -> np.ndarray:
    """降雨強度式 I(t) = a/(t^n+b) [mm/hr] (t=降雨継続時間[分])。
    t=0は定義できないため最小刻み幅程度の正値にクリップする。"""
    t = np.clip(t_min, 1e-6, None)
    return a / (t**n + b)


def generate_central_concentration_hyetograph(
    duration_min: float, dt_min: float, a: float, b: float, n: float,
) -> tuple[np.ndarray, np.ndarray]:
    """降雨強度式I(t)=a/(t^n+b)から, 交互ブロック法により中央集中型ハイエトグラフを作成する。

    各継続時間t_k=k*dt_minでの累積雨量P(t_k)=I(t_k)*t_k/60[mm]の差分ΔP_kを増分雨量とし,
    大きい順に中央から左右交互に配置する(降雨ピークは継続時間の中央付近が標準)。"""
    n_steps = int(round(duration_min / dt_min))
    k = np.arange(1, n_steps + 1)
    t_k = k * dt_min
    I_k = rainfall_intensity_idf(t_k, a, b, n)  # mm/hr
    P_k = I_k * t_k / 60.0  # 累積雨量[mm]
    dP = np.diff(np.concatenate([[0.0], P_k]))  # 各ブロックの増分雨量[mm] (降順)

    rain = np.zeros(n_steps)
    center = n_steps // 2
    positions = [center]
    left, right = center - 1, center + 1
    place_right = True
    while len(positions) < n_steps:
        if place_right and right < n_steps:
            positions.append(right)
            right += 1
        elif not place_right and left >= 0:
            positions.append(left)
            left -= 1
        elif right < n_steps:
            positions.append(right)
            right += 1
        else:
            positions.append(left)
            left -= 1
        place_right = not place_right

    for pos, depth in zip(positions, dP):
        rain[pos] = depth

    t = np.arange(n_steps) * dt_min
    return t, rain


def upsample_rain_to_seconds(t_min: np.ndarray, rain: np.ndarray, dt_min: float) -> tuple[np.ndarray, np.ndarray]:
    """dt_min間隔(既定5分)ブロックの降雨深rainを, 各ブロック内は強度一定とみなして
    1秒刻みに展開する(ブロック内のmm深さをブロック内の秒数で均等分配する)。"""
    sub_steps = int(round(dt_min * UNIT_MIN_TO_SEC))
    rain_1s = np.repeat(rain, sub_steps) / sub_steps
    t_sec_1s = np.arange(len(rain_1s), dtype=float)
    t_min_1s = t_sec_1s / UNIT_MIN_TO_SEC
    return t_min_1s, rain_1s


@dataclass
class ModifiedRRLParams:
    """修正RRL法のパラメータ。

    alpha : (不浸透直接流出, 不浸透間接流出, 浸透直接流出, 浸透間接流出) の寄与率。合計1。
    r_bar : 浸透域の地下浸透成分 R̄ [mm]
    beta  : (beta1, beta2) 凹地貯留量パラメータ [mm]
    """

    alpha: tuple[float, float, float, float] = (0.35, 0.25, 0.25, 0.15)
    r_bar: float = 0.5
    beta: tuple[float, float] = (3.0, 5.0)

    def __post_init__(self) -> None:
        if abs(sum(self.alpha) - 1.0) > 1e-6:
            raise ValueError("alpha の総和は1である必要があります")


def effective_rainfall_modified_rrl(rain: np.ndarray, params: ModifiedRRLParams) -> np.ndarray:
    """修正RRL法による有効降雨量算出。"""
    a1, a2, a3, a4 = params.alpha
    r_bar = params.r_bar
    beta1, beta2 = params.beta

    r_cum = np.cumsum(rain)  # RI(t): 積算降雨量

    def f(r_i: np.ndarray, r: np.ndarray, beta: float) -> np.ndarray:
        # 凹地貯留の増分損失: 累積雨量が容量betaに達するまでは増分雨量rがすべて吸収され,
        # betaを超えたら貯留は満杯のため損失0になる。
        return np.where(r_i > beta, 0.0, r)

    f1 = f(r_cum, rain, beta1)
    f2 = f(r_cum, rain, beta2)

    re = (
        a1 * rain
        + a2 * (rain - f1)
        + a3 * (rain - r_bar)
        + a4 * (rain - r_bar - f2)
    )
    return np.clip(re, 0, None)


def virtual_inflow(effective_rain_mm: np.ndarray, time_area_km2: np.ndarray, dt_min: float) -> np.ndarray:
    """等到達時間域図(Time-Area Curve)との畳込みによる仮想流入量P(t) [m^3/s] の算出。"""
    n = len(effective_rain_mm)
    p_mm_km2 = np.convolve(effective_rain_mm, time_area_km2)[:n]
    # 単位変換: 降雨深1mm x 集水面積1km^2 = 水量1000 m^3
    return p_mm_km2 * 1000.0 / (dt_min * 60.0)


def storage_function_outflow(virtual_inflow_m3s: np.ndarray, K_min: float, dt_min: float) -> np.ndarray:
    """貯留関数法(N=1固定)による流出量Q(t) [m^3/s] の算出。

    dQ/dt = -(1/K) Q + (1/K) P を後退差分近似: Q(t) = ( Q(t-1) + b*P(t) ) / a,  a=1+T/K, b=T/K
    """
    a = 1.0 + dt_min / K_min
    b = dt_min / K_min
    q = np.zeros_like(virtual_inflow_m3s)
    for t in range(1, len(q)):
        q[t] = (q[t - 1] + b * virtual_inflow_m3s[t]) / a
    return q


@dataclass
class CatchmentModel:
    """対象流域の設定。"""

    dt_min: float = 5.0
    time_area_km2: np.ndarray = field(
        default_factory=lambda: np.array([0.45, 0.55, 0.30, 0.20])
    )
    K_min: float = 15.0  # 貯留関数法の貯留定数K[分]

    @property
    def total_area_km2(self) -> float:
        return float(np.sum(self.time_area_km2))


def run_rrl_simulation(rain: np.ndarray, catchment: CatchmentModel,
                        modified_params: ModifiedRRLParams | None = None) -> dict[str, np.ndarray]:
    """降雨データから, 修正RRL法(非線形)+貯留関数法で参考流出量Q_ref(t)を計算する。"""
    params = modified_params or ModifiedRRLParams()
    re = effective_rainfall_modified_rrl(rain, params)
    p = virtual_inflow(re, catchment.time_area_km2, catchment.dt_min)
    q = storage_function_outflow(p, catchment.K_min, catchment.dt_min)
    return {"R": rain, "RE": re, "P": p, "Q": q}


# ----------------------------------------------------------------------
# 2. 降雨 -> 流出量 (MRAC型オンライン適応推定モデル)
#    (appendix_q_model_simulation.py, appendix_q_model_adaptive_control_simulation.py)
# ----------------------------------------------------------------------
@dataclass
class AdaptiveModelConfig:
    p: float                          # P2: 窓積分項の指数 alpha2:=1/P2
    tau_min: float                    # 移動窓の長さ[分]
    theta1_true: float                # 既知の物理定数 theta1:=1/K2 (真の系・推定器で共用)
    theta2_true: float                # 既知の物理定数 theta2:=1/K2 (真の系・推定器で共用)
    p1: float                         # P1: 瞬時減衰項の指数 alpha1:=P1/P2
    catchment_area_km2: float         # 対象流域面積[km2] (I(mm)->I_flow(m3/分)換算に使用)
    theta3_true: float = 1.0          # 真値 theta3=K3 (唯一の未知パラメータ)

    theta3_hat0_ratio: float = 1.5

    gamma3: float = 8.0

    theta3_min: float = 0.05
    theta3_max: float = 5.0

    def __post_init__(self) -> None:
        self.alpha1 = self.p1 / self.p


def _simulate_reference_trajectory(
    I_pos: np.ndarray, dt_min: float, window: int, alpha2: float,
    theta1_true: float, theta2_true: float, theta3_true: float, alpha1: float,
) -> np.ndarray:
    """真値theta1,theta2,theta3で積分した"真の"状態X(t)(=流出量Q(t)^p)。適応則には依存しない。
    PhiI(t) = ∫[t-tau,t] I(s)^alpha2 ds (rain_runoff_delay_replicator_model.pdf 式(19))を
    そのまま使う(推定器側と同一のR(t)を共有するため, 誤差方程式e_dotから降雨項が
    恒等的に消去できる; §1.3参照)。"""
    n = len(I_pos)
    X = np.zeros(n)
    for t in range(1, n):
        start = max(0, t - window)
        PhiI = trapezoid(I_pos[start:t] ** alpha2, dx=dt_min) if t > start else 0.0
        PhiQ = trapezoid(X[start:t] ** alpha2, dx=dt_min) if t > start else 0.0
        target = X[t - 1] + dt_min * (theta1_true * PhiI - theta2_true * PhiQ)
        if target <= 0.0:
            X[t] = 0.0
        else:
            def residual(x: float, target: float = target) -> float:
                return x + dt_min * theta3_true * x**alpha1 - target

            hi = max(target, 1.0)
            while residual(hi) < 0.0:
                hi *= 2.0
            X[t] = brentq(residual, 0.0, hi)
    return X


def _smooth_projection(theta: float, u: float, theta_min: float, theta_max: float,
                        boundary_frac: float = 0.1) -> float:
    """スムーズプロジェクション(標準的なsmooth projection algorithmの簡易版)。

    単純なmin/maxクリップだと, thetaが境界に達した際に更新量uをそのまま加算してから
    事後的に切り詰めるため, e(t)の符号反転が頻発する局面で境界間を高速往復する
    チャタリングが起きうる。境界に近づく(距離がboundary_frac*(theta_max-theta_min)
    未満)につれ境界方向への更新量を連続的に0へ絞ることで, thetaは境界に滑らかに
    漸近するのみとなりチャタリングが起こらない(Lyapunov解析の結論dV/dt<=0を壊さない)。
    """
    span = theta_max - theta_min
    if span <= 0.0:
        return 0.0
    width = boundary_frac * span
    dist_to_max = theta_max - theta
    dist_to_min = theta - theta_min
    if u > 0.0 and dist_to_max < width:
        return u * max(0.0, dist_to_max / width)
    if u < 0.0 and dist_to_min < width:
        return u * max(0.0, dist_to_min / width)
    return u


def _rate_limited_step(du_raw: float, theta_min: float, theta_max: float,
                        max_step_frac: float = 0.2) -> float:
    """1ステップあたりの変化量duを, 可動域(theta_max-theta_min)のmax_step_frac倍
    までに制限する(離散時間Euler近似特有の数値安定化)。

    _smooth_projectionは境界に近づくにつれ更新を絞るが, これは無限小ステップを
    前提とした機構であり, 有限のステップ幅では1ステップの変化量dt_min*uが
    可動域を超えると境界を飛び越えてしまう。変化量自体を可動域の一定割合以下に
    制限することで, thetaが複数ステップかけて滑らかに境界へ近づくよう強制する。
    """
    span = theta_max - theta_min
    if span <= 0.0:
        return 0.0
    max_step = max_step_frac * span
    return max(-max_step, min(max_step, du_raw))


def simulate_adaptive_q_model(
    I: np.ndarray,
    dt_min: float,
    config: AdaptiveModelConfig
) -> dict[str, np.ndarray]:
    """rain_runoff_delay_replicator_model.pdf §1.3のMRAC型適応モデル。

    theta1(=1/K2), theta2(=1/K2)は既知の物理定数として真の系・推定器の両方に
    そのまま使う(オンライン推定するのはtheta3のみ)。

    単位系: 降雨I(t)[mm/区間]をI_flow(t):=I(t)/dt_min*A [m^3/分]の体積流量に直し,
    流出量側の状態XもQ_flow(t):=X(t)^(1/p) [m^3/分]という同じ体積流量として扱う
    ことで, 降雨・流出量が同じ単位系で比較でき, theta1,theta2は面積換算を含まない
    1/K_storageで足りる。戻り値の"Q"は呼び出し側との互換性のため, Q_flow[m^3/分]を
    60で割ってm^3/sへ変換した値である。

    theta3はprojection + rate limit(離散化特有の境界飛び越えを防ぐ数値安定化)を
    かけてEuler積分する(理論式: d(theta3_hat)/dt = -gamma3*e*X^alpha1)。

    注意: dt_minを小さくする(例: 1秒=1/60分)ほどwindow=round(tau_min/dt_min)と
    ステップ数nがともに増え, 各ステップのtrapezoid積分・brentq求根の呼び出し回数が
    増えるため, 計算時間はdt_minに反比例して増大する(5分刻み->1秒刻みで約300倍)。
    """
    if config.tau_min <= 0:
        raise ValueError("tau_min must be positive.")
    if config.p <= 0:
        raise ValueError("p must be positive.")
    if config.theta1_true <= 0:
        raise ValueError("theta1_true must be positive.")
    if config.theta2_true <= 0:
        raise ValueError("theta2_true must be positive.")
    if config.gamma3 <= 0:
        raise ValueError("gamma3 must be positive.")
    if config.catchment_area_km2 <= 0:
        raise ValueError("catchment_area_km2 must be positive.")

    alpha2 = 1.0 / config.p
    alpha1 = config.alpha1

    theta1 = config.theta1_true
    theta2 = config.theta2_true
    theta3_true = config.theta3_true

    # 降雨深[mm]に流域面積[km2]を掛けて体積化(1mm x 1km^2 = 1000m^3)し, 区間長
    # dt_min[分]で割って瞬時流量[m^3/分]にする(積分刻みも分単位のため単位を揃える)。
    I_pos = np.clip(I, 0.0, None)
    I_flow = I_pos / dt_min * config.catchment_area_km2 * UNIT_MM_KM2_TO_M3  # m^3/分

    n = len(I_pos)

    window = int(round(config.tau_min / dt_min))
    if window < 1:
        raise ValueError("tau_min が dt_min に対して小さすぎます(window < 1)。")

    # 真値状態 (PhiI(t) = ∫I_flow(s)^alpha2 ds, 式(19))
    X_true = _simulate_reference_trajectory(
        I_flow, dt_min, window, alpha2, theta1, theta2, theta3_true, alpha1,
    )

    theta3_hat0 = config.theta3_hat0_ratio * theta3_true
    gamma3 = config.gamma3

    X_hat = np.zeros(n)
    theta3_hat = np.zeros(n)
    e = np.zeros(n)

    X_hat[0] = 0.0
    theta3_hat[0] = theta3_hat0
    e[0] = X_true[0] - X_hat[0]

    for t in range(1, n):
        start = max(0, t - window)

        # 移動窓積分 (R(t) = theta1*PhiI(t) - theta2*PhiQ(t), 式(18))。
        # theta1,theta2は既知定数のため, 真の系・推定器の両方に同一のR(t)を使う
        # (PhiQの被積分関数はX_trueを使うseries-parallel構成)。
        if t > start:
            PhiI = trapezoid(I_flow[start:t + 1] ** alpha2, dx=dt_min)
            PhiQ = trapezoid(X_true[start:t + 1] ** alpha2, dx=dt_min)
        else:
            PhiI = 0.0
            PhiQ = 0.0

        R = theta1 * PhiI - theta2 * PhiQ
        th3 = theta3_hat[t - 1]

        # 推定モデル (式(21)): Xhat_dot = -theta3_hat*Xhat^alpha1 + R(t) を
        # 後退Eulerで解く: Xhat[t] + dt*theta3_hat*Xhat[t]^alpha1 = target
        target = X_hat[t - 1] + dt_min * R
        if target <= 0.0:
            X_hat[t] = 0.0
        else:
            def residual(x: float, target: float = target, th3: float = th3) -> float:
                return x + dt_min * th3 * x**alpha1 - target

            hi = max(target, 1.0)
            while residual(hi) < 0.0:
                hi *= 2.0
            X_hat[t] = brentq(residual, 0.0, hi)

        e[t] = X_true[t] - X_hat[t]

        # MRAC回帰信号・適応則 (式(30): theta3_hat_dot = -gamma3*e*X_true^alpha1)
        reg3 = X_true[t] ** alpha1
        theta3_hat_dot = -gamma3 * e[t] * reg3

        # rate limit + スムーズプロジェクションをかけてから加算し, theta3_hatが
        # [theta3_min, theta3_max]の間で境界を飛び越えて往復するチャタリングを防ぐ。
        du3 = _rate_limited_step(dt_min * theta3_hat_dot, config.theta3_min, config.theta3_max)
        du3 = _smooth_projection(theta3_hat[t - 1], du3, config.theta3_min, config.theta3_max)
        theta3_hat[t] = np.clip(theta3_hat[t - 1] + du3, config.theta3_min, config.theta3_max)

    # Q_flow[m^3/分] = X_hat^(1/p) は入力側で既に体積化済みの流量なので, 分->秒の
    # 換算(60で割る)だけでm^3/sになる。
    Q_flow_m3_per_min = X_hat ** (1.0 / config.p)
    Q_m3s = Q_flow_m3_per_min / UNIT_MIN_TO_SEC

    return {
        "I": I_pos,
        "Q": Q_m3s,
        "X_true": X_true,
        "X_hat": X_hat,
        "e": e,
        "theta3_hat": theta3_hat,
        "theta1_true": theta1,
        "theta2_true": theta2,
        "theta3_true": theta3_true,
        "alpha1": alpha1,
    }


def build_adaptive_model_config(cfg: PipelineConfig, catchment_area_km2: float) -> AdaptiveModelConfig:
    """PipelineConfigから, オンライン推定用のAdaptiveModelConfigを組み立てる。

    theta1,theta2は既知の物理定数としてそのまま使う(オンライン推定しない)。
    降雨・流出量を双方ともm^3/分単位の体積流量として扱う(simulate_adaptive_q_model
    参照)ため, 面積換算を個別に含める必要がなく, theta1=theta2=1/K2(=1/K_storage)。"""
    theta1_true = 1.0 / cfg.K_storage
    theta2_true = 1.0 / cfg.K_storage
    return AdaptiveModelConfig(
        p=cfg.p,
        p1=cfg.p1,
        tau_min=cfg.tau_min,
        theta1_true=theta1_true,
        theta2_true=theta2_true,
        catchment_area_km2=catchment_area_km2,
        theta3_true=cfg.adaptive_theta3_true,
        theta3_hat0_ratio=cfg.adaptive_theta3_hat0_ratio,
        gamma3=cfg.adaptive_gamma3,
        theta3_min=cfg.adaptive_theta3_min,
        theta3_max=cfg.adaptive_theta3_max,
    )


def compute_runoff_Q(rain: np.ndarray, dt_min: float, config: AdaptiveModelConfig) -> dict[str, np.ndarray]:
    """降雨I(t)からMRAC型オンライン推定による流出量Q̂(t) [m^3/s] を計算する(戻り値は
    theta3_hat等の時系列を含むsimulate_adaptive_q_modelの結果一式。theta1,theta2は
    既知定数のため定数値のみ)。"""
    return simulate_adaptive_q_model(rain, dt_min=dt_min, config=config)


def resample_to_1s(t_min: np.ndarray, Q_m3s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """dt_min刻みのQ(t)を1秒刻みに線形補間し, 「1秒ごとに更新される流入量」を模擬する。"""
    t_sec = t_min * UNIT_MIN_TO_SEC
    t_sec_1s = np.arange(t_sec[0], t_sec[-1] + 1.0, 1.0)
    Q_1s = np.interp(t_sec_1s, t_sec, Q_m3s)
    return t_sec_1s, Q_1s


# ----------------------------------------------------------------------
# 3. 降雨-流出量ピーク遅れ時間の推定 (次回tau_min設定の判断材料)
#    (replicator_sinc_simulation_delay_embedding.py, rain_runoff_delay_replicator.py)
# ----------------------------------------------------------------------
def sinc_kernel_matrix(v: np.ndarray) -> np.ndarray:
    """ペアワイズsincカーネル行列: K[i,j] = sinc(v_i - v_j)。"""
    diff = v[:, None] - v[None, :]
    return np.sinc(diff)


def _build_delay_embedding_rhs(I_norm, Q_norm, delays_min: np.ndarray):
    """I_norm(t)(降雨), Q_norm(t)(流出量, ピーク値で正規化した補間関数)から,
    降雨側を遅延埋め込みし, 候補遅延k*tauだけ過去の降雨I(t-k*tau)が現在の流出量
    Q(t)にどれだけsinc類似度で一致するかをタップkの適応度とするfitness/replicator_rhs
    を組み立てて返す。

    delays_min: 候補遅延タップの配列(呼び出し側でK_min[分]以下を除外して渡す)。"""

    def xt_vec(t):
        return I_norm(t - delays_min)

    def fitness(t):
        return sinc_kernel_matrix(np.append(xt_vec(t), Q_norm(t)))[:-1, -1]

    def replicator_rhs(t, omega):
        f = fitness(t)
        payoff = omega @ f
        return omega * (f - payoff)

    return xt_vec, fitness, replicator_rhs


def estimate_peak_delay(t_min: np.ndarray, rain: np.ndarray, Q_5min: np.ndarray,
                         label: str, out_path: str, K_min: float,
                         cfg: PipelineConfig) -> tuple[float, float]:
    """降雨rain(t_min)と流出量Q_5min(t_min)から, レプリケータ動力学によるピーク遅れ時間を推定する。

    候補遅延タップはcfg.delay_taps_minのうちK_min[分]以下を除外して使う(推定される
    遅れ時間tlagは常にK_minを超え, 次回tau設定の判断材料tau_suggest=2*(tlag-K_min)は
    必ず0以上になる。採否は出力グラフ・数値を見て人が判断し, 自動反映はしない)。
    K_min以下しかタップが残らない場合は, K_minを超える範囲まで自動延長する。

    最小の候補タップ遅延delays_min[0]分より前は, どのタップの降雨信号も過去データの
    範囲外(0埋め)で情報を持たないため, レプリケータの計算・表示から除外する。

    戻り値: (推定遅れ時間tlag[分], 次回tauの目安tau_suggest[分])。"""
    delays_min_full = cfg.delay_taps_min
    delays_min = delays_min_full[delays_min_full > K_min]
    if len(delays_min) == 0:
        first_tap_above = cfg.delay_tap_min * (int(K_min // cfg.delay_tap_min) + 1)
        delays_min = first_tap_above + cfg.delay_tap_min * np.arange(cfg.n_taps)
        print(f"  [情報] 候補遅延タップ(0〜{delays_min_full[-1]:.0f}分)がすべてK_min({K_min:.1f}分)"
              f"以下だったため, {delays_min[0]:.0f}〜{delays_min[-1]:.0f}分の範囲に自動延長しました。")
    n_taps_eff = len(delays_min)

    I_max = float(np.max(rain))
    Q_max = float(np.max(Q_5min))
    I_interp = interp1d(t_min, rain / I_max, bounds_error=False, fill_value=0.0)
    Q_interp = interp1d(t_min, Q_5min / Q_max, bounds_error=False, fill_value=0.0)

    xt_vec, fitness, replicator_rhs = _build_delay_embedding_rhs(I_interp, Q_interp, delays_min)

    T_total = float(t_min[-1])
    # t=delays_min[0]より前は情報を持たないため, その時刻から計算・表示を開始する。
    T_start = float(delays_min[0])
    omega0 = np.ones(n_taps_eff) / n_taps_eff
    t_eval = np.arange(T_start, T_total + cfg.dt_eval, cfg.dt_eval)
    sol = solve_ivp(
        replicator_rhs, t_span=(T_start, T_total), y0=omega0, t_eval=t_eval,
        method=cfg.ivp_method, max_step=cfg.ivp_max_step,
        rtol=cfg.ivp_rtol, atol=cfg.ivp_atol,
    )

    omega_final = sol.y[:, -1]
    i_star = int(np.argmax(omega_final))
    delay_est_min = float(delays_min[i_star])
    if i_star == n_taps_eff - 1:
        print(f"  [警告] 最大遅延タップ({delays_min[-1]:.0f}分)に張り付いています。"
              f"n_tapsを増やして再実行してください。")

    peak_I_t = float(t_min[np.argmax(rain)])
    peak_Q_t = float(t_min[np.argmax(Q_5min)])

    tau_suggest_min = 2.0 * (delay_est_min - K_min)

    print(f"  レプリケータ推定 遅れ時間: {delay_est_min:.0f}分  "
          f"(タップ{i_star}/{n_taps_eff - 1}, omega={omega_final[i_star]:.3f}, "
          f"K_min={K_min:.0f}分以下を{len(delays_min_full) - n_taps_eff}タップ除外済み)")
    print(f"  (参考) 降雨ピーク時刻={peak_I_t:.0f}分, 流出量ピーク時刻={peak_Q_t:.0f}分")
    # tau_suggestは全イベント処理後にmain()側でまとめて提案として表示する。

    # 表示用: シフトなしの生の降雨I_norm(t)・流出量Q_norm(t)。
    I0 = np.array([I_interp(t) for t in t_eval])
    Qn0 = np.array([Q_interp(t) for t in t_eval])
    fitness_star = np.array([fitness(t)[i_star] for t in t_eval])

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)

    axes[0].plot(t_eval, I0, label="降雨 I_norm(t)(ピーク値で正規化)")
    axes[0].plot(t_eval, Qn0, label="流出量 Q_norm(t)(ピーク値で正規化)")
    axes[0].axvline(peak_I_t, color="tab:blue", linestyle=":", label=f"降雨ピーク t={peak_I_t:.0f}分")
    axes[0].axvline(peak_Q_t, color="tab:orange", linestyle=":", label=f"流出量ピーク t={peak_Q_t:.0f}分")
    axes[0].legend(loc="upper right", fontsize=8)
    axes[0].set_ylabel("正規化信号")
    axes[0].set_title(f"{label}: 降雨・流出量波形(ともに自身のピーク値で正規化)")

    for i in range(n_taps_eff):
        show_label = i in (0, i_star, n_taps_eff - 1)
        axes[1].plot(
            sol.t, sol.y[i],
            label=rf"$\omega_{{{i}}}$ (遅延{delays_min[i]:.0f}分)" if show_label else None,
            linewidth=2.0 if i == i_star else 0.8,
            alpha=1.0 if i == i_star else 0.5,
        )
    axes[1].legend(loc="upper left", fontsize=8)
    axes[1].set_ylabel(r"$\omega_i(t)$")
    axes[1].set_title(f"レプリケータ重み (N={n_taps_eff}, K_min={K_min:.0f}分以下除外, "
                       f"{cfg.delay_tap_min:g}分刻み) -> 推定遅れ時間={delay_est_min:.0f}分")

    axes[2].plot(t_eval, fitness_star)
    axes[2].set_ylabel(rf"fitness$_{{{i_star}}}(t)$")
    axes[2].set_xlabel("時刻 [分]")
    axes[2].set_title(f"推定タップ(遅延{delay_est_min:.0f}分)のsinc適応度 = sinc(I(t-{delay_est_min:.0f}分) - Q(t))")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"プロットを保存しました: {out_path}")

    return delay_est_min, tau_suggest_min


# ----------------------------------------------------------------------
# 4. ポンプ場物理モデル (resources_07_00073_pump_scheduling_simulation.py)
# ----------------------------------------------------------------------
@dataclass
class PumpParams:
    """可変速水中ポンプ1台の性能モデル。"""

    Nmax: float = 30.0
    Nmin: float = 15.0
    Q_bep_Nmax: float = 0.20
    H_bep_Nmax: float = 12.0
    H_shutoff_ratio: float = 1.25
    eta_bep_Nmax: float = 0.70
    eta_bep_drop_at_Nmin: float = 0.15
    rel_eff_curvature: float = 0.7
    rel_eff_floor: float = 0.15

    def __post_init__(self) -> None:
        # 揚程曲線 H_N(Q) = c0h*N^2 + c2h*Q^2 を「締切揚程(Q=0)」と「BEP点」の2点で較正する。
        H_shutoff_Nmax = self.H_shutoff_ratio * self.H_bep_Nmax
        self._c0h = H_shutoff_Nmax / self.Nmax**2
        self._c2h = (self.H_bep_Nmax - H_shutoff_Nmax) / self.Q_bep_Nmax**2
        self._x_bep = self.Q_bep_Nmax / self.Nmax

    def head(self, Q, N):
        return self._c0h * N**2 + self._c2h * Q**2

    def eta_bep(self, N):
        slope = self.eta_bep_drop_at_Nmin * self.eta_bep_Nmax / (self.Nmax - self.Nmin)
        return self.eta_bep_Nmax - slope * (self.Nmax - N)

    def relative_efficiency(self, Q, N):
        x = Q / N
        e = 1.0 - self.rel_eff_curvature * ((x - self._x_bep) / self._x_bep) ** 2
        return np.clip(e, self.rel_eff_floor, 1.0)

    def efficiency(self, Q, N):
        return self.relative_efficiency(Q, N) * self.eta_bep(N)

    def power(self, Q, N):
        H = np.maximum(self.head(Q, N), 0.0)
        eta = np.maximum(self.efficiency(Q, N), 1.0e-6)
        return GAMMA_WATER * H * Q / eta


@dataclass
class PlantScenario:
    """ウェットウェル+配管シナリオの無次元パラメータ alpha, beta で指定する。"""

    S: float = 10.0
    beta: float = 0.5
    alpha: float = 1.5
    Sh_max: float = 10.0
    Qad_max_ratio: float = 3.0


@dataclass
class PlantDerived:
    """PlantScenario + PumpParams から導出される, シミュレーションで直接使う量。"""

    S: float
    H0: float
    K: float
    Qin_max: float
    Hw_min: float
    Hw_max: float
    Qad_max: float


def derive_plant(pump: PumpParams, scenario: PlantScenario) -> PlantDerived:
    H0 = scenario.beta * pump.H_bep_Nmax
    K = (pump.H_bep_Nmax - H0) / pump.Q_bep_Nmax**2
    Qin_max = pump.Q_bep_Nmax / scenario.alpha
    W = 900.0 * pump.Q_bep_Nmax / scenario.Sh_max
    Hw_max = W / scenario.S
    Qad_max = scenario.Qad_max_ratio * pump.Q_bep_Nmax
    return PlantDerived(S=scenario.S, H0=H0, K=K, Qin_max=Qin_max, Hw_min=0.0,
                         Hw_max=Hw_max, Qad_max=Qad_max)


def energy_kwh(P: np.ndarray, dt_s: float) -> float:
    return float(np.sum(P) * dt_s / 3.6e6)


def reference_energy_kwh(Qin: np.ndarray, plant: PlantDerived, dt_s: float) -> float:
    """損失を無視した理論上の最小ポンプ電力量Eref。"""
    integrand = GAMMA_WATER * Qin * (plant.H0 + plant.K * Qin**2)
    return float(trapezoid(integrand, dx=dt_s) / 3.6e6)


def build_pump_station(cfg: PipelineConfig) -> tuple[PumpParams, PlantDerived]:
    """元論文の小型ポンプ1台想定のPumpParams/PlantScenarioを, 適応推定モデル側の
    流出量Q(t)の規模(峰値数十〜二百数十m^3/s)に合わせて大型雨水ポンプ場の規模に
    合わせ直す(単位は[m^3/s]のまま)。

    Hw_max=900*Q_bep_Nmax/(Sh_max*S)をcfg.plant_Hw_max_targetにしたいので,
    S=900*Q_bep_Nmax/(Sh_max*Hw_max_target)から逆算する。
    """
    pump = PumpParams(
        Nmax=cfg.pump_Nmax, Nmin=cfg.pump_Nmin,
        Q_bep_Nmax=cfg.pump_Q_bep_Nmax,
        H_bep_Nmax=cfg.pump_H_bep_Nmax,
        H_shutoff_ratio=cfg.pump_H_shutoff_ratio,
        eta_bep_Nmax=cfg.pump_eta_bep_Nmax,
        eta_bep_drop_at_Nmin=cfg.pump_eta_bep_drop_at_Nmin,
        rel_eff_curvature=cfg.pump_rel_eff_curvature,
        rel_eff_floor=cfg.pump_rel_eff_floor,
    )
    S = 900.0 * cfg.pump_Q_bep_Nmax / (cfg.plant_Sh_max * cfg.plant_Hw_max_target)  # m^2
    scenario = PlantScenario(S=S, beta=cfg.plant_beta, alpha=1.0,
                              Sh_max=cfg.plant_Sh_max, Qad_max_ratio=cfg.plant_Qad_max_ratio)
    plant = derive_plant(pump, scenario)
    return pump, plant


# ----------------------------------------------------------------------
# 5. 流出量 -> ポンプ運転台数 (複数台ポンプの段階起動)
#
# 同一機種のポンプn_pumps台が, 共通のウェットウェル・吐出管で並列運転する構成。
# 水位に応じて1台ずつ段階的に起動/停止し(ヒステリシス制御), 合計消費電力を求める。
# ----------------------------------------------------------------------
def solve_operating_point_multi(Hw: float, N: float, n_online: int,
                                 pump: PumpParams, plant: PlantDerived) -> float:
    """n_online台(>=1)を並列運転したときの合計吐出量Qtotalを求める。n_online<=0なら0.0。"""
    if n_online <= 0:
        return 0.0
    denom = pump._c2h / n_online**2 - plant.K
    numer = plant.H0 - Hw - pump._c0h * N**2
    Q2 = numer / denom
    if Q2 <= 0.0 or np.isnan(Q2):
        return 0.0
    Q_total = float(np.sqrt(Q2))
    return float(min(Q_total, plant.Qad_max))


def build_stage_thresholds(plant: PlantDerived, n_pumps: int) -> np.ndarray:
    """[Hw_min, Hw_max]をn_pumps等分した(n_pumps+1)個の水位閾値を返す。
    k台目(k=1..n_pumps)は水位がedges[k]以上で起動し, edges[k-1]以下で停止する。"""
    return np.linspace(plant.Hw_min, plant.Hw_max, n_pumps + 1)


def cs_control_step_multi(Hw_prev: float, m_prev: int, Qin_i: float,
                           pump: PumpParams, plant: PlantDerived,
                           edges: np.ndarray, dt_s: float):
    """複数台版の1秒ステップ。前の秒の水位Hw_prevで台数mを決め, 流入と吐出を同時反映する。
    戻り値: (Hw, m, Q_total, N, P_total)。"""
    n_pumps = len(edges) - 1
    m = m_prev
    while m < n_pumps and Hw_prev >= edges[m + 1]:
        m += 1
    while m > 0 and Hw_prev <= edges[m - 1]:
        m -= 1

    if m > 0:
        Q_total = solve_operating_point_multi(Hw_prev, pump.Nmax, m, pump, plant)
        N = pump.Nmax
    else:
        Q_total = 0.0
        N = 0.0

    Hw_new = Hw_prev + (Qin_i - Q_total) / plant.S * dt_s
    Hw = float(np.clip(Hw_new, plant.Hw_min, plant.Hw_max))

    if m > 0:
        Q_each = Q_total / m
        P_total = m * float(pump.power(Q_each, N))
    else:
        P_total = 0.0
    return Hw, m, Q_total, N, P_total


def run_realtime_stepwise_multi(t_sec: np.ndarray, Qin: np.ndarray, pump: PumpParams,
                                 plant: PlantDerived, n_pumps: int,
                                 dt_s: float = 1.0) -> dict[str, np.ndarray]:
    """cs_control_step_multiを1秒ごとに逐次呼び出し, n_pumps台構成の全区間をシミュレートする。"""
    edges = build_stage_thresholds(plant, n_pumps)
    n = len(Qin)
    Hw = np.zeros(n)
    Q = np.zeros(n)
    N = np.zeros(n)
    P = np.zeros(n)
    m_online = np.zeros(n, dtype=int)
    # 水位はウェットウェルが空の状態(Hw_min)から開始する(半分溜まった状態から
    # 始めると, 流入量が小さい場合に架空の初期貯留水を排出しようとして開始直後に
    # 実態と無関係な吐出スパイクが発生するため)。
    Hw[0] = plant.Hw_min
    m = 0
    for i in range(1, n):
        Hw[i], m, Q[i], N[i], P[i] = cs_control_step_multi(Hw[i - 1], m, Qin[i], pump, plant, edges, dt_s)
        m_online[i] = m
    return {"t_sec": t_sec, "Hw": Hw, "Q": Q, "N": N, "P": P, "m_online": m_online, "edges": edges}


# ----------------------------------------------------------------------
# 6. 可視化
# ----------------------------------------------------------------------
def evaluate_mrac_stability_conditions(
    cfg: PipelineConfig, adaptive_config: AdaptiveModelConfig,
) -> list[tuple[str, bool, str]]:
    """このMRAC型適応推定モデルが理論上安定(漸近安定AS, 追従誤差e(t)->0)であるための
    条件を, 現在のPipelineConfigの値で評価する(theta1,theta2は既知定数のため,
    オンライン推定するtheta3に関する条件のみが残る)。

    条件: P1>0,P2>0(モデルの適切性), P1>P2(名目系AS優位), P1<=1(摂動下でもAS),
    theta3_min>0(Barbalatの補題による収束にprojection下限の正値が必要),
    gamma3>0(適応ゲイン符号), theta1>0/theta2>0(既知定数が正)。

    戻り値: [(条件の説明, 満たすか, 現在値の詳細), ...]。"""
    P1, P2 = cfg.p1, cfg.p
    return [
        (
            "P1>0, P2>0 (適切性: Q^P1, X:=Q^P2の定義)",
            P1 > 0 and P2 > 0,
            f"P1={P1:g}, P2={P2:g}",
        ),
        (
            "P1>P2 (式(P1-P2-condition): nu=P1/P2-1>0, 名目系AS優位)",
            P1 > P2,
            f"P1={P1:g}, P2={P2:g}",
        ),
        (
            "P1<=1 (式(rho-nu-condition): 摂動R(t)が仮定2を満たしAS)",
            P1 <= 1.0,
            f"P1={P1:g}",
        ),
        (
            "theta3_min>0 (projection下限が正: Barbalatの補題によるe(t)->0の結論に必要)",
            cfg.adaptive_theta3_min > 0,
            f"theta3_min={cfg.adaptive_theta3_min:g}",
        ),
        (
            "gamma3>0 (適応ゲイン正)",
            cfg.adaptive_gamma3 > 0,
            f"gamma3={cfg.adaptive_gamma3:g}",
        ),
        (
            "theta1>0, theta2>0 (既知定数theta1=theta2:=1/K2が正)",
            adaptive_config.theta1_true > 0 and adaptive_config.theta2_true > 0,
            f"theta1={adaptive_config.theta1_true:g}, theta2={adaptive_config.theta2_true:g}",
        ),
    ]


def plot_theta_estimates(t_min: np.ndarray, res: dict[str, np.ndarray], cfg: PipelineConfig,
                          adaptive_config: AdaptiveModelConfig, title: str, out_path: str) -> None:
    """オンライン推定(MRAC型適応則)によるthetâ3(t)の収束の様子を表示する(theta1,theta2は
    既知定数のため推定せず, theta3のみ表示)。あわせてthetâ3の上下限(projection範囲)と,
    理論上の安定条件(evaluate_mrac_stability_conditions参照)の充足状況を左上に表示する。"""
    conditions = evaluate_mrac_stability_conditions(cfg, adaptive_config)
    all_ok = all(ok for _, ok, _ in conditions)

    n_cond = len(conditions)
    top_margin = 0.10 + 0.028 * n_cond
    fig, ax_th3 = plt.subplots(1, 1, figsize=(10, 4.5 + 2.2 * top_margin))

    ax_th3.axhspan(cfg.adaptive_theta3_min, cfg.adaptive_theta3_max, color="tab:green", alpha=0.08,
                    label=f"projection範囲[{cfg.adaptive_theta3_min:.3g}, {cfg.adaptive_theta3_max:.3g}]")
    ax_th3.plot(t_min, res["theta3_hat"], color="tab:green", linewidth=1.6, label=r"$\hat{\theta}_3(t)$")
    ax_th3.axhline(res["theta3_true"], color="black", linestyle="--", linewidth=1.2, label="真値 theta3=K3")
    ax_th3.set_title(f"{title}: オンライン推定値thetâ3の収束(theta1,theta2は既知定数)")
    ax_th3.set_xlabel("時刻 [分]")
    ax_th3.set_ylabel(r"$\hat{\theta}_3(t)$")
    ax_th3.grid(alpha=0.3)
    ax_th3.legend(fontsize=8)

    status_lines = ["rain_runoff_delay_replicator_model.pdfに基づく安定条件:"]
    for label_text, ok, detail in conditions:
        status_lines.append(f"  [{'OK' if ok else 'NG'}] {label_text}  ({detail})")
    box_color = "#e6f4ea" if all_ok else "#fdecea"
    edge_color = "tab:green" if all_ok else "tab:red"
    fig.text(0.01, 0.995, "\n".join(status_lines), ha="left", va="top", fontsize=8,
              bbox=dict(boxstyle="round", facecolor=box_color, edgecolor=edge_color, alpha=0.92))

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 1.0 - top_margin))
    fig.savefig(out_path, dpi=150)
    print(f"プロットを保存しました: {out_path}")


def plot_state_X_1s(t_min_1s: np.ndarray, res_1s: dict[str, np.ndarray], title: str, out_path: str) -> None:
    """1秒刻みでシミュレートした状態X(t)(=Q(t)^p)の真値X_trueと推定値X_hatを表示する。"""
    fig, ax = plt.subplots(figsize=(12, 5))
    t_h = t_min_1s / 60.0
    ax.plot(t_h, res_1s["X_true"], color="black", linewidth=1.0, label="真値 X_true(t)")
    ax.plot(t_h, res_1s["X_hat"], color="tab:green", linewidth=1.0, label="推定値 X_hat(t)")
    ax.set_xlabel("時刻 [時]")
    ax.set_ylabel(r"状態 X(t) = Q(t)^p")
    ax.set_title(f"{title}: 状態X(t)の時系列(1秒刻みシミュレーション)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"プロットを保存しました: {out_path}")


def plot_pump_schedule(t_min: np.ndarray, Q_5min: np.ndarray, Q_ref: np.ndarray,
                        t_sec: np.ndarray, Qin_1s: np.ndarray,
                        sim: dict[str, np.ndarray], plant: PlantDerived, n_pumps: int,
                        title: str, out_path: str) -> None:
    """流量・稼働台数m(t)・水位・消費電力の4段グラフ。"""
    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
    ax_q, ax_m, ax_hw, ax_p = axes
    t_h = t_sec / 3600.0

    ax_q.plot(t_min / 60.0, Q_5min, "o", color="gray", markersize=3,
              label="流出量 Q̂(t) (オンライン推定値, 5分刻み)")
    ax_q.plot(t_min / 60.0, Q_ref, ":", color="dimgray", linewidth=1.6,
              label="Q_ref(t) (参考: 修正RRL法+貯留関数法, 5分刻み)")
    ax_q.plot(t_h, Qin_1s, "-", color="tab:blue", linewidth=1.0, label="流入量Qin(t) (1秒刻みシミュレーション)")
    ax_q.plot(t_h, sim["Q"], "-", color="tab:green", linewidth=1.0, label="ポンプ合計吐出量Qtotal(t)")
    ax_q.axhline(plant.Qad_max, color="tab:red", linestyle="--", linewidth=1.0, label="ポンプ最大吐出量Qad_max(閾値)")
    ax_q.set_ylabel("流量 [m3/s]")
    ax_q.set_title(f"{title}  (同一機種ポンプ{n_pumps}台の段階起動)")
    ax_q.grid(alpha=0.3)
    ax_q.legend(fontsize=8)

    ax_m.step(t_h, sim["m_online"], where="post", color="tab:brown", linewidth=1.4, label="稼働台数 m(t)")
    ax_m.set_ylabel(f"稼働台数 [台] (0〜{n_pumps})")
    ax_m.set_yticks(range(n_pumps + 1))
    ax_m.grid(alpha=0.3)
    ax_m.legend(fontsize=8)

    ax_hw.plot(t_h, sim["Hw"], color="tab:green", linewidth=1.2, label="ウェットウェル水位 Hw(t)")
    for k, edge in enumerate(sim["edges"]):
        style = "--" if k in (0, n_pumps) else ":"
        ax_hw.axhline(edge, color="gray", linestyle=style, linewidth=0.8)
    ax_hw.set_ylabel("水位 [m]\n(点線=各台の起動/停止閾値)")
    ax_hw.grid(alpha=0.3)
    ax_hw.legend(fontsize=8)

    ax_p.plot(t_h, sim["P"] / 1000.0, color="tab:purple", linewidth=1.2, label="合計消費電力 P_total(t)")
    ax_p.set_xlabel("時刻 [時]")
    ax_p.set_ylabel("消費電力 [kW]")
    ax_p.grid(alpha=0.3)
    ax_p.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"プロットを保存しました: {out_path}")


def plot_rain_runoff_pump(t_min: np.ndarray, rain: np.ndarray, Q_5min: np.ndarray,
                           t_sec: np.ndarray, Q_pump: np.ndarray, title: str, out_path: str) -> None:
    """降雨量I(t)・流出量Q̂(t)・ポンプ合計吐出量Qtotal(t)を1枚に重ねたグラフ。
    降雨量は単位・桁が流量と異なるため, 右軸反転の棒グラフで表示する。"""
    fig, ax = plt.subplots(figsize=(12, 6))
    t_h_min = t_min / 60.0
    t_h_sec = t_sec / 3600.0
    dt_min = t_min[1] - t_min[0]

    ax.plot(t_h_min, Q_5min, color="tab:red", linewidth=1.8, label="流出量 Q̂(t) (オンライン推定値)")
    ax.plot(t_h_sec, Q_pump, color="tab:green", linewidth=1.2, label="ポンプ合計吐出量 Qtotal(t)")
    ax.set_xlabel("時刻 [時]")
    ax.set_ylabel("流量 [m3/s]")
    ax.set_title(f"{title}: 降雨量・流出量・ポンプ吐出量")
    ax.grid(alpha=0.3)

    ax_r = ax.twinx()
    ax_r.bar(t_h_min, rain, width=(dt_min / 60.0) * 0.8, color="tab:blue", alpha=0.4, label="降雨量 I(t)")
    ax_r.set_ylabel("降雨量 [mm]")
    ax_r.invert_yaxis()

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax_r.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"プロットを保存しました: {out_path}")


# ----------------------------------------------------------------------
# 7. メイン処理: 降雨 -> 流出量 -> ポンプ運転台数
#
# 5分刻みの通常パイプライン(RRL参考流量・遅れ時間推定)に加え, 状態X(t)を1秒刻みで
# 見るため, simulate_adaptive_q_modelをdt_min=1/60(=1秒)で再実行する。window・
# ステップ数が5分刻み比で約300倍になるため, 計算時間も大幅に増える(コメント参照)。
# ----------------------------------------------------------------------
def main() -> None:
    cfg = PipelineConfig()

    catchment = CatchmentModel(dt_min=cfg.catchment_dt_min, time_area_km2=cfg.catchment_time_area_km2,
                                K_min=cfg.catchment_K_min)
    adaptive_config = build_adaptive_model_config(cfg, catchment.total_area_km2)
    pump, plant = build_pump_station(cfg)
    rrl_params = ModifiedRRLParams()

    print(
        f"対象流域面積: A={catchment.total_area_km2:.2f}km2\n"
        f"ポンプ場規模: Q_bep_Nmax={pump.Q_bep_Nmax:.1f}m3/s, H_bep_Nmax={pump.H_bep_Nmax:.1f}m, "
        f"Nmax={pump.Nmax:.0f}rpm/Nmin={pump.Nmin:.0f}rpm\n"
        f"ウェットウェル: S={plant.S:.1f}m2, Hw_min={plant.Hw_min:.2f}m, Hw_max={plant.Hw_max:.2f}m, "
        f"Qad_max={plant.Qad_max:.1f}m3/s\n"
    )

    events = {
        label: generate_central_concentration_hyetograph(
            duration_min=cfg.rain_duration_min, dt_min=catchment.dt_min, a=a, b=b, n=n,
        )
        for label, (a, b, n) in cfg.rain_idf_params.items()
    }

    # tau_minは実行中固定値として使う。各イベントのtau_suggestはここに集め,
    # 全イベント処理後に次回への提案としてまとめて表示する(main()末尾参照)。
    tau_suggestions: list[tuple[str, float, float]] = []

    for label, (t_min, rain) in events.items():
        print(f"=== {label} ===")

        # --- 1) 降雨 -> 流出量 (MRAC型オンライン適応推定, 5分刻み) ---
        res = compute_runoff_Q(rain, dt_min=catchment.dt_min, config=adaptive_config)
        Q_5min = res["Q"]
        print(f"  流出量Q_hat(t): 峰値={Q_5min.max():.2f}m3/s, 平均={Q_5min.mean():.2f}m3/s "
              f"({len(Q_5min)}点, dt={catchment.dt_min}分)")

        # 比較用: 修正RRL法+貯留関数法(オンライン推定を使わない, 別モデルによる参考流出量)
        result_ref = run_rrl_simulation(rain, catchment, modified_params=rrl_params)
        Q_ref = result_ref["Q"]
        print(f"  参考Q_ref(t)(修正RRL法+貯留関数法): 峰値={Q_ref.max():.2f}m3/s, 平均={Q_ref.mean():.2f}m3/s")

        theta_out_path = Path(__file__).resolve().parent / f"{OUTPUT_PREFIX}_theta_{label[:6]}.png".replace(" ", "")
        plot_theta_estimates(t_min, res, cfg, adaptive_config, label, str(theta_out_path))

        # --- 2) 降雨-流出量ピーク遅れ時間の推定 (次回tau_min設定の判断材料, 5分刻み) ---
        # tau_suggest=2*(tlag-K)のKは適応モデル側のK_storage(catchment.K_minとは別物)。
        delay_out_path = Path(__file__).resolve().parent / f"{OUTPUT_PREFIX}_delay_{label[:6]}.png".replace(" ", "")
        delay_est_min, tau_suggest_min = estimate_peak_delay(
            t_min, rain, Q_5min, label, str(delay_out_path), cfg.K_storage, cfg,
        )
        tau_suggestions.append((label, delay_est_min, tau_suggest_min))

        # --- 3) 状態X(t)を1秒刻みで見るため, シミュレーション自体を1秒刻みで再実行する ---
        # 降雨を5分ブロック内一定強度として1秒刻みに展開し, dt_min=1/60でsimulate_adaptive_
        # q_modelを再実行する(window/ステップ数とも約300倍になり計算時間も同程度増える)。
        t_min_1s, rain_1s = upsample_rain_to_seconds(t_min, rain, catchment.dt_min)
        print(f"  [1秒刻み] {len(rain_1s)}ステップのシミュレーションを開始します(時間がかかります)...")
        t0 = time.perf_counter()
        res_1s = compute_runoff_Q(rain_1s, dt_min=1.0 / UNIT_MIN_TO_SEC, config=adaptive_config)
        elapsed = time.perf_counter() - t0
        print(f"  [1秒刻み] 完了: {elapsed:.1f}秒")

        x_out_path = Path(__file__).resolve().parent / f"{OUTPUT_PREFIX}_stateX_{label[:6]}.png".replace(" ", "")
        plot_state_X_1s(t_min_1s, res_1s, label, str(x_out_path))

        # --- 4) 流出量 -> ポンプ運転台数 (複数台段階起動) ---
        # ポンプ入力Qin_1sは, 補間ではなく1秒刻みシミュレーション結果res_1s["Q"]を直接使う。
        t_sec = t_min_1s * UNIT_MIN_TO_SEC
        Qin_1s = res_1s["Q"]
        sim = run_realtime_stepwise_multi(t_sec, Qin_1s, pump, plant, n_pumps=cfg.n_pumps, dt_s=cfg.dt_s)

        assert not np.any(np.isnan(sim["Hw"])) and not np.any(np.isnan(sim["P"]))
        assert np.all((sim["Hw"] >= plant.Hw_min - 1e-9) & (sim["Hw"] <= plant.Hw_max + 1e-9))
        assert np.all((sim["m_online"] >= 0) & (sim["m_online"] <= cfg.n_pumps))
        assert np.all(sim["P"] >= 0.0)

        E_pump = energy_kwh(sim["P"], dt_s=cfg.dt_s)
        E_ref = reference_energy_kwh(Qin_1s, plant, dt_s=cfg.dt_s)
        n_starts = int(np.sum(np.diff(sim["m_online"]) > 0))
        n_saturated = int(np.sum(Qin_1s > plant.Qad_max))
        print(f"  [ポンプ運転台数({cfg.n_pumps}台構成)] 起動イベント数: {n_starts}回, "
              f"最大同時稼働台数: {sim['m_online'].max()}台")
        print(f"  流入量がポンプ容量Qad_maxを超えた秒数: {n_saturated}/{len(Qin_1s)}"
              f"  (超えた分はポンプが追従できず水位上昇・オーバーフローの恐れ)")
        print(f"  合計消費電力量: {E_pump:.1f}kWh (理論最小Eref: {E_ref:.1f}kWh)\n")

        pump_out_path = Path(__file__).resolve().parent / f"{OUTPUT_PREFIX}_schedule_{label[:6]}.png".replace(" ", "")
        plot_pump_schedule(t_min, Q_5min, Q_ref, t_sec, Qin_1s, sim, plant, cfg.n_pumps, label, str(pump_out_path))

        overlay_out_path = Path(__file__).resolve().parent / f"{OUTPUT_PREFIX}_overlay_{label[:6]}.png".replace(" ", "")
        plot_rain_runoff_pump(t_min, rain, Q_5min, t_sec, sim["Q"], label, str(overlay_out_path))

    # 全イベント処理後にtau_min設定の判断材料としてまとめて表示する(自動反映はしない)。
    print(f"=== 次回の降雨イベントに向けたtau_min設定の提案(現在の固定値: tau_min={cfg.tau_min:.1f}分) ===")
    for label, tlag, tau_suggest in tau_suggestions:
        print(f"  {label}: tlag={tlag:.0f}分 -> tau_suggest=2*(tlag-K_storage)={tau_suggest:.1f}分")
    tau_suggest_mean = float(np.mean([t for _, _, t in tau_suggestions]))
    print(f"  [参考] 全イベント平均のtau_suggest: {tau_suggest_mean:.1f}分  "
          f"(採用するかはこの数値・上記グラフを見て判断すること。自動反映はしません)")

    plt.show()


if __name__ == "__main__":
    main()
