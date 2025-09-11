# -*- coding: utf-8 -*-
# 地図クリック or 緯度経度入力 で地点指定 / 結果保持 / AMD_Tools4 対応
import streamlit as st
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as md

import AMD_Tools4 as AMD                       # AMD_Tools4（USER/PASSWORDS 設定が必要）
from AMD_PaddyWaterTemp import WaterTemp       # WaterTemp は AMD_Tools4 参照版を使用

import folium
from streamlit_folium import st_folium

st.markdown("### 🌾 水田水温取得アプリ（信大作成）")
st.title("🌾 水田水温取得アプリ")

# -------------------------------
# セッション初期化
# -------------------------------
if "clicked_lat" not in st.session_state:
    st.session_state.clicked_lat = 35.7915215
if "clicked_lon" not in st.session_state:
    st.session_state.clicked_lon = 137.9473239
if "result_df" not in st.session_state:
    st.session_state.result_df = None
if "result_meta" not in st.session_state:
    st.session_state.result_meta = None

# -------------------------------
# 入力フォーム（地図とは独立・送信ボタンで確定）
# -------------------------------
with st.sidebar:
    st.header("📍 位置の指定方法")
    mode = st.radio("地点指定の方法", ["地図で選ぶ", "緯度経度を直接入力"], horizontal=False)

    if mode == "緯度経度を直接入力":
        # 直接入力 UI（直近のクリック結果を初期表示に）
        lat_input = st.number_input("緯度 (lat)", value=float(st.session_state.clicked_lat), format="%.6f")
        lon_input = st.number_input("経度 (lon)", value=float(st.session_state.clicked_lon), format="%.6f")
    else:
        lat_input = None
        lon_input = None

    st.header("📅 計算条件")
    with st.form(key="controls"):
        # デフォルト日付を 2025-05-01 ～ 2025-09-30 に設定
        start_date = st.date_input("開始日", datetime(2025, 5, 1))
        end_date   = st.date_input("終了日", datetime(2025, 9, 30))

        water_depth = st.number_input("水深 (mm)", value=50, step=10, help="全期間一定。日別入力の拡張も可能です。")
        dvs_bgn = st.number_input("DVS開始（移植=0）", value=0.0, step=0.1, min_value=0.0, max_value=2.0)
        dvs_end = st.number_input("DVS終了（成熟=2）", value=2.0, step=0.1, min_value=0.0, max_value=2.0)
        submitted = st.form_submit_button("計算実行")

    # 明示的クリア
    if st.button("🧹 結果をクリア"):
        st.session_state.result_df = None
        st.session_state.result_meta = None
        st.success("結果をクリアしました。")

# -------------------------------
# 地図（クリックで座標更新）
# ※ 「地図で選ぶ」モードのときだけ表示
# -------------------------------
if mode == "地図で選ぶ":
    st.subheader("🗺️ 地図をクリックして地点を選択")
    m = folium.Map(
        location=[st.session_state.clicked_lat, st.session_state.clicked_lon],
        zoom_start=12, tiles="OpenStreetMap", control_scale=True
    )
    folium.Marker(
        [st.session_state.clicked_lat, st.session_state.clicked_lon],
        tooltip="選択中の地点",
        icon=folium.Icon(color="green")
    ).add_to(m)
    map_state = st_folium(m, width=None, height=420)
    if map_state and map_state.get("last_clicked") is not None:
        st.session_state.clicked_lat = float(map_state["last_clicked"]["lat"])
        st.session_state.clicked_lon = float(map_state["last_clicked"]["lng"])

# 現在の選択座標（モードに応じて決定）
if mode == "緯度経度を直接入力":
    current_lat = float(lat_input)
    current_lon = float(lon_input)
else:
    current_lat = float(st.session_state.clicked_lat)
    current_lon = float(st.session_state.clicked_lon)

col_a, col_b = st.columns(2)
with col_a:
    st.metric("緯度 (lat)", f"{current_lat:.6f}")
with col_b:
    st.metric("経度 (lon)", f"{current_lon:.6f}")

# -------------------------------
# 計算関数（キャッシュ）
# -------------------------------
@st.cache_data(show_spinner="計算中...", ttl=0)
def run_water_temp(timedomain, lalodomain, dvs, Dw):
    # 実データ（観測/予報）
    Tw_mea, Tw_max, Tw_min, tim, lat0, lon0 = WaterTemp(timedomain, lalodomain, dvs, Dw)
    # 平年値（cli=True、風速 ws が必須 → 固定2 m/sを使用）
    Tw_mea_N, Tw_max_N, Tw_min_N, _, _, _   = WaterTemp(timedomain, lalodomain, dvs, Dw, ws=2, cli=True)

    # 3D → 1D（単一点）
    Tw_mea, Tw_max, Tw_min       = Tw_mea[:,0,0], Tw_max[:,0,0], Tw_min[:,0,0]
    Tw_mea_N, Tw_max_N, Tw_min_N = Tw_mea_N[:,0,0], Tw_max_N[:,0,0], Tw_min_N[:,0,0]

    df = pd.DataFrame({
        "Date": tim,
        "Tw_mea": Tw_mea, "Tw_max": Tw_max, "Tw_min": Tw_min,
        "Tw_mea_N": Tw_mea_N, "Tw_max_N": Tw_max_N, "Tw_min_N": Tw_min_N
    })
    return df

# -------------------------------
# フォーム送信時のみ計算
# -------------------------------
if submitted:
    if end_date <= start_date:
        st.error("終了日は開始日より後にしてください。")
    else:
        timedomain = [start_date.strftime("%Y-%m-%d"), end_date.strftime("%Y-%m-%d")]
        # WaterTemp の仕様：単一点は [lat, lat, lon, lon] で渡す
        lalodomain = [current_lat, current_lat, current_lon, current_lon]  # 単一点指定（仕様に準拠）  # noqa
        # days と DVS・Dw 長さを一致
        days = len(AMD.timedom(timedomain))
        if days <= 0:
            st.error("期間に日数がありません。開始日・終了日を見直してください。")
        else:
            dvs = np.linspace(float(dvs_bgn), float(dvs_end), days, endpoint=False, dtype=float).tolist()
            Dw  = [float(water_depth)] * days
            try:
                df = run_water_temp(timedomain, lalodomain, dvs, Dw)
                st.session_state.result_df = df
                st.session_state.result_meta = {
                    "lat": current_lat, "lon": current_lon,
                    "start": timedomain[0], "end": timedomain[1],
                    "mode": mode
                }
                st.success("計算が完了しました。結果は保持されます。")
            except Exception as e:
                st.error(f"❌ エラー: {e}")
                st.info("※ AMD_Tools4.py の USER/PASSWORDS 設定や、AMD_PaddyWaterTemp.py の配置（AMD_Tools4 参照）をご確認ください。")
                st.session_state.result_df = None
                st.session_state.result_meta = None

# -------------------------------
# 結果表示（保持）
# -------------------------------
if st.session_state.result_df is not None:
    df = st.session_state.result_df
    meta = st.session_state.result_meta or {}
    st.subheader("📊 計算結果（時系列）")
    st.dataframe(df, use_container_width=True)

    csv = df.to_csv(index=False).encode("utf-8-sig")
    st.download_button("📥 CSVをダウンロード", csv, "water_temp_result.csv", "text/csv")

    # プロット（Tw_mea vs Tw_mea_N）
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(df["Date"], df["Tw_mea"], label="日平均水温", color="k")
    ax.plot(df["Date"], df["Tw_mea_N"], label="平年値（日平均）", color="k", linewidth=0.8)
    ax.fill_between(df["Date"], df["Tw_mea"], df["Tw_mea_N"], where=df["Tw_mea"]>df["Tw_mea_N"], alpha=0.5)
    ax.fill_between(df["Date"], df["Tw_mea_N"], df["Tw_mea"], where=df["Tw_mea"]<df["Tw_mea_N"], alpha=0.5)
    ax.set_xlabel("Date"); ax.set_ylabel("Water Temp (°C)")
    if meta:
        ax.set_title(f"水田水温（lat={meta.get('lat',0):.5f}, lon={meta.get('lon',0):.5f} / {meta.get('start','')}–{meta.get('end','')}）")
    ax.xaxis.set_major_locator(md.DayLocator(bymonthday=[1]))
    ax.xaxis.set_major_formatter(md.DateFormatter('%m/%d'))
    ax.legend()
    st.pyplot(fig)
else:
    st.info("サイドバーで条件を設定し、「計算実行」を押してください。結果は画面に保持されます。")


