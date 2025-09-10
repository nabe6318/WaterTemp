# -*- coding: utf-8 -*-
# Streamlit app: Paddy field water temperature (AMD_Tools4 + AMD_PaddyWaterTemp)
# - DVSは移植(0)→成熟(2)を等間隔で進行させる簡易モデル（オプションで閾値を出穂=1固定）
# - 水深 Dw は期間中一定値（mm）
# - 実測/再解析ベースと「平年値」(cli=True) の両系列を計算
# - CSV ダウンロードとグラフ表示、地図で地点確認
#
# 必要モジュール:
#   AMD_Tools4, AMD_PaddyWaterTemp, streamlit, numpy, pandas, matplotlib, folium, streamlit_folium
#
# 作成: 2025-09-11 (for Prof. Watanabe)

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from datetime import datetime, date
import AMD_Tools4 as amd
from AMD_PaddyWaterTemp import WaterTemp
import folium
from streamlit_folium import st_folium

st.set_page_config(page_title="水田水温（AMD_Tools4）", layout="wide")
st.title("🌾 水田水温推定アプリ（AMD_Tools4 + AMD_PaddyWaterTemp）")

st.markdown("""
**機能**  
- 期間・地点・水深を指定して、水田の **日平均・最高・最低水温** を推定  
- **平年値（cli=True）** との比較を同時計算（風速 `ws` は任意指定）  
- **CSV ダウンロード** と比較グラフ表示  
""")

# ---------------------------
# プリセット地点（必要に応じて編集）
# ---------------------------
locations = {
    "北海道岩見沢市（サンプル）": (43.2463, 141.6946),
    "AFC 野辺山ステーション": (35.945, 138.476),
    "伊那キャンパス": (35.827, 137.961),
}

with st.sidebar:
    st.header("▶ 入力設定")

    # 1) 地点選択 or 手入力
    mode = st.radio("地点の指定方法", ["プリセットから選択", "緯度経度を手入力"], index=0)

    if mode == "プリセットから選択":
        loc_name = st.selectbox("地点を選択", list(locations.keys()), index=0)
        lat, lon = locations[loc_name]
    else:
        lat = st.number_input("緯度 (lat)", value=43.2463, format="%.6f")
        lon = st.number_input("経度 (lon)", value=141.6946, format="%.6f")
        loc_name = f"Lat={lat:.4f}, Lon={lon:.4f}"

    # 2) 期間指定（デフォルト: 2017-06-01〜2017-09-28）
    st.markdown("---")
    st.write("**計算期間（日別）**")
    d_beg = st.date_input("開始日", value=date(2017, 6, 1))
    d_end = st.date_input("終了日", value=date(2017, 9, 28))
    if d_end < d_beg:
        st.error("終了日は開始日以降にしてください。")

    # 3) 発育ステージ（DVS）と水深
    st.markdown("---")
    st.write("**生育パラメータ**")
    dvs_bgn = st.number_input("開始DVS（移植=0）", value=0.0, step=0.1)
    dvs_end = st.number_input("終了DVS（成熟=2）", value=2.0, step=0.1)
    enforce_heading = st.checkbox("出穂（DVS=1）を期間中央に固定（等間隔補正）", value=True)

    Dw_mm = st.number_input("水深 Dw [mm]（期間一定）", value=50, step=5, min_value=1)
    ws_cli = st.number_input("平年値計算用の風速 ws [m/s]", value=2.0, step=0.5, min_value=0.1)

    # 4) 解析オプション
    st.markdown("---")
    st.write("**解析オプション**")
    make_cumsum = st.checkbox("積算（累積）表示に切り替える", value=False)

# 地図で位置確認
with st.expander("📍 地図で地点を確認（Folium）", expanded=False):
    m = folium.Map(location=[lat, lon], zoom_start=10)
    folium.Marker([lat, lon], tooltip=loc_name).add_to(m)
    st_folium(m, width=600, height=400)

# ---------------------------
# 計算ボタン
# ---------------------------
calc = st.button("水田水温を計算")

if calc:
    try:
        # --- 時間軸（日別） ---
        timedomain = [str(d_beg), str(d_end)]
        # AMD_Tools4 の日別タイムドメイン（datetime の配列を受けとる想定）
        td = amd.timedom(timedomain)  # e.g., [datetime(YYYY,MM,DD), ...]
        days = len(td)
        if days <= 0:
            st.error("期間内の日数がゼロです。日付を見直してください。")
            st.stop()

        # --- DVS を等間隔で用意 ---
        # 例: np.linspace(dvs_bgn, dvs_end, days, endpoint=False) として日数ぶんの等間隔列
        dvs = np.linspace(dvs_bgn, dvs_end, days, endpoint=False)

        # 出穂(DVS=1)を期間中央に合わせたい場合の簡易補正（線形スケーリング）
        if enforce_heading and dvs_end > dvs_bgn:
            # 期間中央インデックス
            mid_idx = days // 2
            # 現在の中点DVS
            cur_mid = dvs[mid_idx]
            # 全体をシフト・スケールして mid が 1 になるよう補正
            # スケールは元の幅を維持、まずシフトで mid を 1 へ
            shift = 1.0 - cur_mid
            dvs = dvs + shift
            # 始端・終端が外れすぎる場合は上下でクリップ
            dvs = np.clip(dvs, 0.0, 2.0)

        # --- 水深 Dw（一定値） ---
        Dw = [float(Dw_mm)] * days

        # --- LatLonDomain 形式（latmin, latmax, lonmin, lonmax）---
        lalodomain = [lat, lat, lon, lon]

        # --- 水温計算（実系列／平年値系列）---
        # 実系列（ws: デフォルトは1kmメッシュから）
        Tw_mea, Tw_max, Tw_min, tim, lat_a, lon_a = WaterTemp(
            timedomain, lalodomain, dvs, Dw
        )
        # 平年値系列（cli=True、風速は ws_cli を与える）
        Tw_mea_N, Tw_max_N, Tw_min_N, timN, _, _ = WaterTemp(
            timedomain, lalodomain, dvs, Dw, ws=float(ws_cli), cli=True
        )

        # --- 3D -> 1D に変換 ---
        def squeeze3(x):
            x = np.asarray(x)
            if x.ndim == 3:
                return x[:, 0, 0]
            return x

        Tw_mea = squeeze3(Tw_mea)
        Tw_max = squeeze3(Tw_max)
        Tw_min = squeeze3(Tw_min)
        Tw_mea_N = squeeze3(Tw_mea_N)
        Tw_max_N = squeeze3(Tw_max_N)
        Tw_min_N = squeeze3(Tw_min_N)

        # 時刻を pandas.DatetimeIndex へ
        tim = pd.to_datetime(tim)

        # --- 積算モード（任意） ---
        M_mea = Tw_mea.copy()
        M_meaN = Tw_mea_N.copy()
        if make_cumsum:
            M_mea = np.cumsum(M_mea)
            M_meaN = np.cumsum(M_meaN)

        # --- CSV 出力用 DataFrame ---
        df = pd.DataFrame({
            "Date": tim.date,
            "Tw_mea": Tw_mea,
            "Tw_max": Tw_max,
            "Tw_min": Tw_min,
            "Tw_mea_N": Tw_mea_N,
            "Tw_max_N": Tw_max_N,
            "Tw_min_N": Tw_min_N,
            "DVS": dvs,
            "Dw_mm": Dw_mm,
            "Lat": lat,
            "Lon": lon
        })

        st.subheader("📄 計算結果（表）")
        st.dataframe(df, use_container_width=True)

        # ダウンロードボタン
        csv = df.to_csv(index=False).encode("utf-8-sig")
        st.download_button(
            label="CSV をダウンロード",
            data=csv,
            file_name="paddy_water_temperature_AMD4.csv",
            mime="text/csv"
        )

        # --- グラフ（実系列 vs 平年値, 日平均 or 積算） ---
        st.subheader("📈 比較グラフ")
        fig, ax = plt.subplots(figsize=(12, 4))

        # 影塗り（実系列が平年値より高い/低い日）
        ax.fill_between(tim, M_mea, M_meaN, where=(M_mea > M_meaN), alpha=0.5)
        ax.fill_between(tim, M_meaN, M_mea, where=(M_mea < M_meaN), alpha=0.5)

        ax.plot(tim, M_mea, 'k', label="実系列（日平均）" if not make_cumsum else "実系列（累積）")
        ax.plot(tim, M_meaN, 'k', linewidth=0.6, label="平年値（日平均）" if not make_cumsum else "平年値（累積）")

        ax.set_xlabel("日付")
        ax.set_ylabel("水温 [℃]" if not make_cumsum else "累積水温 [℃・日]")
        ax.set_title(f"{loc_name}  N{lat:.4f}, E{lon:.4f}")
        ax.xaxis.set_major_formatter(md.DateFormatter('%m/%d'))
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.legend(loc="best")
        plt.tight_layout()
        st.pyplot(fig)

        with st.expander("詳細パラメータ・注記"):
            st.markdown(f"""
- 計算日数: **{days} 日**
- DVS: {dvs_bgn} → {dvs_end}（等間隔）, 出穂固定: **{enforce_heading}**
- 水深: **{Dw_mm} mm**（期間一定）
- 平年値の風速 ws: **{ws_cli} m/s**
- `AMD_Tools4.timedom()` により日単位の時系列を生成  
- `AMD_PaddyWaterTemp.WaterTemp()` を用いて **Tw_mea/Tw_max/Tw_min** を推定  
- 平年値は `cli=True` で取得（風速は平年値がないため **ws** 指定）
""")

    except Exception as e:
        st.error(f"計算エラー: {e}")
        st.stop()
