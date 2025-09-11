# -*- coding: utf-8 -*-
"""
PaddyWaterTempの関数「WaterTemp」の説明

WaterTemp関数
    WaterTemp関数は、ユーザが任意で指定した開始日と終了日、緯度・経度範囲、イネの発育ステージ、最大葉面積指数、水深、アルベド、風速の入力値から
            気象データをもとに水田の日平均・最高・最低水温を出力します。気象データはAMD_Tools3.pyを利用して1kmメッシュ農業気象データサーバから取得します。
            風速は固定値を直接指定することもできます。
            
書式：
    WaterTemp(timedomain, lalodomain, dvs, Dw, laim=4.6, alb=0.1, ws=None, cli=False)

引数（必須）：
    timedomain：計算する日付を範囲で、['2008-05-05', '2008-05-06']のような文字列の2要素リストで与える。
                特定の日を計算するときは、二か所に同じ日付を与える。
    lalodomain：計算するデータの緯度と経度の範囲で、[36.0, 40.0, 130.0, 135.0]のように南端緯度、北端緯度、西端経度、東端経度の順で指定する。
                特定地点のデータを計算するときは、緯度と経度にそれぞれ同じ値を与える。
    dvs：期間中の水稲の発育ステージを[0.1, 0.2, 0.3]のような日数分のリストで与える。発育ステージの定義は0:移植,1:出穂,2:成熟。
    Dw：水深を[50, 40, 30]のような日数分のリストで与える。単位はmm。
    
引数（必要に応じ指定）：
    laim：最大葉面積指数を数値で与える。省略すると4.6が設定される。日々の葉面積指数LAIはdvsとlaimから計算される。dvs=0とした場合はLAI=0。単位はm2/m2。
    alb：アルベドを数値で与える。省略すると0.1が設定される。開始日から終了日まで一定。単位は無次元。
    ws：風速を数値で与える。開始日から終了日まで一定。省略すると1kmメッシュ農業気象データサーバから取得する。単位はm/s
    cli：平年値を取得するときにcli=True（または1）として指定する。省略した場合は観測値や予報値を使用する。ただし風速には平年値がないのでwsで設定する。
    
戻り値：
    第1戻り値：日平均水温（浮動小数の三次元（日付×緯度×経度）配列）。単位は℃。
    第2戻り値：日最高水温（浮動小数の三次元（日付×緯度×経度）配列）。単位は℃。
    第3戻り値：日最低水温（浮動小数の三次元（日付×緯度×経度）配列）。単位は℃。
    第4戻り値：計算した日付の並び（日時オブジェクトの一次元配列）。
    第5戻り値：計算したメッシュの中心緯度の並び（浮動小数の一次元配列）。
    第6戻り値：計算したメッシュの中心経度の並び（浮動小数の一次元配列）。
    
使用例：緯度35度、経度135度の地点の2017年1月1日～2017年1月3日の水田水温を計算する場合。
    2日間の発育ステージは0.1、0.2、0.3、水深は50mm、40mm、30mmと変化させ、風速は1kmメッシュ農業気象サーバのデータを使用する。
    from PaddyWaterTemp import WaterTemp
    timedomain = ['2017-01-01', '2017-01-03']  #計算期間の設定
    lalodomain  = [35, 35, 135, 135]  #計算領域の設定
    dvs = [0.1, 0.2, 0.3] 	#期間中の発育ステージ(期間分のリストであること)
    Dw = [50, 40, 30] 	#期間中の水深(期間分のリストであること)
    Tmea, Tmax, Tmin, tim, lat, lon = WaterTemp(timedomain, lalodomain, dvs, Dw)
 
Atsushi MARUYAMA, 2018.01.16
"""

from __future__ import print_function

from datetime import datetime as dt
from datetime import timedelta as td
from argparse import ArgumentParser
import numpy as np
from math import pi,log,exp,atan,asin,acos,log10,sqrt,sin,cos,tan
from cmath import exp as cexp
# ▼▼▼ ここを AMD_Tools4 に変更 ▼▼▼
import AMD_Tools4 as AMD
# ▲▲▲ ここを AMD_Tools4 に変更 ▲▲▲

GSR_FACTOR = 10000/864
LAI = 1e-3
h   = 2e-3
DVS = 1e-3
MESHSIZE = 5068800

#####################################################################################
#       Physiological constant

RCd = 287.0   # Gas constant of the dry air (J kg^-1 K^-2)
RCw = 461.5   # Gas constant of the water vapor (J kg^-1 K^-2)
cpd = 1005.0  # Spesific heat of the dry air (J kg^-1 K^-2)
cpw = 1854.0  # Spesific heat of the water vapor (J kg^-1 K^-2)
P0 = 1013.2   # Standard atmospheric pressure (hPa)

Gravitation = 9.8       #Graviational acceleration (m s^-2)
Karman = 0.4            #Karman constant (ND)
STB = 5.67e-8           #Stefan-Boltzmann constant(W m^-2 K^-4)
I00 = 1365.0            #Solar constant(W m^-2)

##############################################################################
# Transcendental function

def sec(x): return 1/cos(x)

##############################################################################
# Empirical and theoretical functions for the air and water vapor

def cal_Patm(Elevation):
    '''Empirical: Patm(hPa), Elevation(m)'''
    return P0 - 0.1093 * Elevation

def esat(T):
    '''esat(hPa), T(C)'''
    A = 7.5
    B = 237.3
    if T < -200: T = -200
    return 6.1078 * 10 ** (A * T / (B + T))

def cal_e(T, RH):
    '''e(hPa), T(C), RH(%)'''
    return esat(T) * RH / 100

def cal_q(e, Patm):
    '''q(kg kg^-1), e(hPa), Patm(hPa)'''
    return 0.622 * (e / Patm) / (1 - 0.378 * (e / Patm))

def cal_e_from_q(q, Patm):
    '''e(hPa), q(kg kg^-1), Patm(hPa)'''
    return Patm * q / (0.622 + 0.378 * q)

def qsat(T, Patm):
    '''qsat(kg kg^-1), T(C), Patm(hPa)'''
    return cal_q(esat(T), Patm)

def cal_RC(q):
    return  RCd * (1 - q) + RCw * q

def cal_cp(q):
    return cpd * (1 - q) + cpw * q

def cal_rho(T, e, Patm):
    '''rho(kg m^-3), T(C), e(hPa), Patm(hPa)'''
    return 1.293 * (273.15 / (273.15 + T)) * (Patm / P0) * (1 - 0.378 * (e / Patm))

def cal_Latent(T):
    '''Latent Heat(J kg^-1), T(C)'''
    return 2.5 * 10 ** 6 - 2400 * T

def cal_theta(T, RC, cp, Patm):
    '''potential temp(C), Patm(hPa) 
       RC(J kg^-1 K^-2), cp(J kg^-1 K^-2)'''
    return (T + 273.15) * (1000 / Patm) ** (RC / cp) - 273.15

def BlackRad(T):
    '''BlackRad(W m^-2), T(C)'''
    return STB * (T + 273.15) ** 4

###############################################################################
#   Estimation of radiation

def cal_Tdew(e):
    '''e(hPa), Tdew(C)'''
    A = 7.5
    B = 237.3
    if e == 0.0: return -B
    return B * log10(e / 6.1078) / (A - log10(e / 6.1078))
    
def cal_SolarAngle(DOY, JST, Latitude, Longitude):
    DOYrad = DOY * (2 * pi / 365)
    Sekii = asin(0.398 * sin(4.871 + DOYrad + 0.033 * sin(DOYrad)))
    Akashi_Longitude = 135 * (2 * pi / 360)                     # Standard longitude in Japan
    NanchuTime = 12 - (Longitude - Akashi_Longitude) * 24 / (2 * pi)
    Jikaku_from_nanchu = (JST - NanchuTime) * (2 * pi) / 24
    return asin(sin(Latitude) * sin(Sekii) +\
                cos(Latitude) * cos(Sekii) * cos(Jikaku_from_nanchu))    # Elevation angle

def cal_S0day(DOY, Latitude):
    DOYrad = DOY * (2 * pi / 365)
    Sekii = asin(0.398 * sin(4.871 + DOYrad + 0.033 * sin(DOYrad)))
    DistanceEffect = (1.00011 + 0.034221 * cos(DOYrad) + 0.00128 * sin(DOYrad) +\
                      0.000719 * cos(DOYrad * 2) + 0.000077 * sin(DOYrad * 2))
    Jikaku = acos(-tan(Latitude) * tan(Sekii))
    S0day = I00 / pi * DistanceEffect *\
            (Jikaku * sin(Latitude) * sin(Sekii) +\
             cos(Latitude) * cos(Sekii) * sin(Jikaku))
    N0 = Jikaku * 24 / pi
    return Sekii, S0day, N0

def cal_Sday(S0day, N, N0):
    A = 0.244
    B = 0.511
    C = 0.118
    if N == 0:
        Sday = S0day * C
    else:
        Sday = S0day * (A + B * N / N0)
    return Sday

def cal_Sfday(Latitude, Sekii, S0day, ea, Patm):
    #print("Latitude, Sekii, S0day, ea, Patm",Latitude, Sekii, S0day, ea, Patm)
    Bdust = 0.03
    refwide = 0.15
    #print("LAT",Latitude,Sekii,sec(Latitude - Sekii))
    log10w = 0.0312 * cal_Tdew(ea) - 0.0963
    C1 = max(0.21 - 0.2 * Bdust, 0.15)
    F1 = 0.056 + 0.16 * sqrt(Bdust)
    j1 = (0.066 + 0.34 * sqrt(Bdust)) * (refwide - 0.15)
    mnoon = Patm / P0 * sec(Latitude - Sekii)
    k3 = 1.402 - 0.06 * log10(Bdust + 0.02) - 0.1 * sqrt(sec(Latitude - Sekii) - 0.91)
    mday = k3 * mnoon
    i3 = 0.014 * (mday + 7 + 2 * log10w) * log10w
    Sfday = S0day * (C1 + 0.7 * 10 ** (-mday * F1)) * (1 - i3) * (1 + j1)
    #print("log10w,C1,F1,j1,mnoon,k3,mday,i3,Sfday",log10w,C1,F1,j1,mnoon,k3,mday,i3,Sfday)
    return Sfday

def cal_Lday(Ta, ea, Ceffect):
    log10wtop = 0.0315 * cal_Tdew(ea) - 0.1836
    Lfday = (0.74 + 0.19 * log10wtop + 0.07 * log10wtop ** 2) * BlackRad(Ta)
    Lday = BlackRad(Ta) * (1 - (1 - Lfday / BlackRad(Ta)) * Ceffect)
    return Lday

def CeffectA(N, N0):
    A = N / N0
    return 0.826 * A ** 3 - 1.234 * A ** 2 + 1.135 * A + 0.298


def CeffectB(Sday, Sfday):
    B = Sday / Sfday
    return 0.03 * B ** 3 - 0.3 * B ** 2 + 1.25 * B - 0.04

#'------------------------------------------------------------------------------
def cal_Lday_fromN(DOY, Latitude, Ta, ea, N):
    Sekii, S0day, N0 = cal_S0day(DOY, Latitude)
    Sday = cal_Sday(S0day, N, N0)
    Ceffect = CeffectA(N, N0)
    Lday = cal_Lday(Ta, ea, Ceffect)
    return (N0, Sday, Lday)


def cal_Lday_fromS(DOY, Latitude, Ta, ea, Patm, Sday):
    #print("DOY, Latitude, Ta, ea, Patm, Sday",DOY, Latitude, Ta, ea, Patm, Sday)
    Sekii, S0day, N0 = cal_S0day(DOY, Latitude)
    #print("Sekii, S0day, N0",Sekii, S0day, N0)
    Sfday = cal_Sfday(Latitude, Sekii, S0day, ea, Patm)
    #print("Sfday",Sfday)
    Ceffect = CeffectB(Sday, Sfday)
    #print("Ceffect",Ceffect)
    Lday = cal_Lday(Ta, ea, Ceffect)
    #print("Lday",Lday)
    return (N0, Lday)


###############################################################################
#   Diurnal variation in air tempetaure and solar radiation

def cal_Ta(Tave, Tmax, Tmin, Time24):
    B01 = -1 * (Tmax - Tmin) / 2.09
    B02 = -0.2 * B01
        
    return Tave + B01 * cos(1 * (2 * pi / 24) * Time24 - pi / 4) \
        + B02 * cos(2 * (2 * pi / 24) * Time24 - pi / 4)

    
def cal_Sd(Sday, Time24):
    A01 = -1.503
    A02 = 0.584
    A03 = -0.058
    A04 = -0.023
    
    return Sday * (1 + A01 * cos(1 * (2 * pi / 24) * Time24) \
                       + A02 * cos(2 * (2 * pi / 24) * Time24) \
                       + A03 * cos(3 * (2 * pi / 24) * Time24) \
                       + A04 * cos(4 * (2 * pi / 24) * Time24))



###############################################################################
#   Diurnal variation in wind speed (2012/12 add)

def cal_Ua(Uave, Time24):
    C01 = -0.32
    C02 = 0.093
    
    return Uave * (1 + C01 * cos(1 * (2 * pi / 24) * Time24 - pi / 4) \
                       + C02 * cos(2 * (2 * pi / 24) * Time24 - pi / 4))
 

###############################################################################
# Double Source Model   2006/03/29 by Atsushi Maruyama
#                       2011/07/13 Modification of stability effects
###############################################################################

def DSM(Patm, Ta, ea, Ua, Sd, Ld, h, LAI, DVS, za, Tg0, TgM, Dw, refg):

    #------------------------------------------------------------------------------
    # Parameters of the model (Maruyama et al. [2017] Submitted)
    h_geo_s = 0.025     # Geometrical roughness of the soil surface (m)
    Dw_f = 50           # Water depth where soil surface is submerged (mm)
    c_free_w = 0.0024   # Heat exchange coefficient by free convection on water surface (ND)
    c_free_s = 0.0036   # Heat exchange coefficient by free convection on soil surface (ND)
    cdf = 0.2           # Drag coefficient of the leaf (ND)
    chf = 0.06          # Heat transfer coefficient of the leaf (ND)
    alpha = 0.4         # Absorptivity of the leaf for short-wave radiation (ND)
    # Site dependent parameters
    #refg = 0.1          # Albedo of the ground (ND)
    #------------------------------------------------------------------------------

    TimeStep = 3600
    TimeCycle = TimeStep * 24
    omega = 2 * pi / TimeCycle

    qa = cal_q(ea, Patm)
    RC = cal_RC(qa)
    cp = cal_cp(qa)
    rho = cal_rho(Ta, ea, Patm)
    Latent = cal_Latent(Ta)
    theta = cal_theta(Ta, RC, cp, Patm)

    F = cal_F(DVS)

    ref, tauS, tauL = cal_RadPara(LAI, F, alpha, refg)

    Sg = (1 - refg) * tauS * Sd
    Sc = (1 - ref) * Sd - Sg
    Scmean = Sc / float(LAI)
    Ucmean = cal_Ucmean(za, h, LAI, cdf, Ua)

    gsf = cal_gsf(DVS, Scmean)
    cef = gsf * chf / (gsf + chf * Ucmean)
    if cef < 0.0001: cef = 0.0001

    # ------ Calculation of rouhness and c_free from water depth ------

    h_geo = h_geo_s * max(1 - Dw / Dw_f, 0)
    c_free = c_free_w + (c_free_s - c_free_w) * (h_geo / h_geo_s)

    h_geo = max(h_geo, 0.00001)
    ustar_geo = 0.03    # Friction velocity for the ground(m/s)
    z0s, zTs, zqs = cal_Roughness(za, h_geo, ustar_geo)

    #------------------------------------------------------------------

    Tgmax = Ta + 30
    Tgmin = Ta - 30
    Tcmax = Ta + 30
    Tcmin = Ta - 30
    
    count1 = 0

    d, z0, zT, zq, z0t, zTt, zqt = cal_AeroPara(h, LAI, z0s, zTs, zqs, cdf, chf, cef)

    #------------------------------------------------------------------
    # Calculation of bulk tranfer coefficients under neutral condition
    #------------------------------------------------------------------
    CM = cal_CX(za, d, z0, z0, 0, 0)
    CH = cal_CX(za, d, z0, zT, 0, 0)
    CE = cal_CX(za, d, z0, zq, 0, 0)

    CMg = cal_CX(za, d, z0, z0t, 0, 0)
    CHg = cal_CX(za, d, z0, zTt, 0, 0)
    CEg = cal_CX(za, d, z0, zqt, 0, 0)
    
    CMs = cal_CX(za, 0, z0s, z0s, 0, 0)
    CHs = cal_CX(za, 0, z0s, zTs, 0, 0)
    CEs = cal_CX(za, 0, z0s, zqs, 0, 0)

    fmg = CMg / float(CMs)
    fhg = CHg / float(CHs)
    feg = CEg / float(CEs)


    fhc = (CH - CHg) / float(CH)
    fec = (CE - CEg) / float(CE)

    psiM = 0
    psiH = 0
    psiE = 0
    psiMg = 0
    psiHg = 0
    psiEg = 0

    #----------------------------------------------------------------------
    # Calculation of bulk tranfer coefficients under non-neutral condition
    #----------------------------------------------------------------------
    while True:
        CM = cal_CX(za, d, z0, z0, psiM, psiM)
        CH = cal_CX(za, d, z0, zT, psiM, psiH)
        CE = cal_CX(za, d, z0, zq, psiM, psiE)

        CMs = cal_CX(za, 0, z0s, z0s, psiMg, psiMg)
        CHs = cal_CX(za, 0, z0s, zTs, psiMg, psiHg)
        CEs = cal_CX(za, 0, z0s, zqs, psiMg, psiEg)

        CMg = fmg * CMs
        CHg = fhg * CHs
        CEg = feg * CEs

        CHc = fhc * CH
        CEc = fec * CE

        count1 = count1 + 1
        if count1 >= 99: break      # In case when energy budget is not convergent

        Tg = (Tgmax + Tgmin) / 2.0  # Bisection method
        Tc = (Tcmax + Tcmin) / 2.0  # Bisection method

        #********* Flux calculation ****************************************************************

        CHUg = max(CHg * Ua, c_free * max((Tg - Ta), 0) ** (1 / 3.0))  # Heat exchange by free convection (2011/04)
        CEUg = max(CEg * Ua, c_free * max((Tg - Ta), 0) ** (1 / 3.0))  # Heat exchange by free convection (2011/04)

        HHg = cp * rho * CHUg * (Tg - Ta)
        LEg = Latent * rho * CEUg * (qsat(Tg, Patm) - qa)

        HHc = cp * rho * CHc * Ua * (Tc - Ta)
        LEc = Latent * rho * CEc * Ua * (qsat(Tc, Patm) - qa)

        HH = HHg + HHc
        LE = LEg + LEc

        Rng = Sg + tauL * Ld + (1 - tauL) * BlackRad(Tc) - BlackRad(Tg)
        Rnc = Sc + (1 - tauL) * (Ld + BlackRad(Tg) - 2 * BlackRad(Tc))
        Rn = Rng + Rnc

        G = sqrt(omega * 1000 * 3430 * 1.21 / 2.0) * ((Tg - Tg0) / float(omega * TimeStep) + Tg - TgM) # Force-restore model (for soil layer)
        G = G + (4180 * (Tg - Tg0) / float(TimeStep) * Dw)         # Heat flux into water layer

        #*******************************************************************************************

        if Rng - (HHg + LEg + G) < 0:
            Tgmax = (Tg + Tgmax) / 2.0
        else:
            Tgmin = (Tg + Tgmin) / 2.0

        if Rnc - (HHc + LEc) < 0:
            Tcmax = (Tc + Tcmax) / 2.0
        else:
            Tcmin = (Tc + Tcmin) / 2.0

        psiM, psiH, psiE = cal_psi(theta, CM, Ua, HH, cp, rho, za, d, z0, zT, zq)
        psiMg, psiHg, psiEg = cal_psi(theta, CMg, Ua, HHg, cp, rho, za, 0, z0s, zTs, zqs)

        if abs(Rng - (HHg + LEg + G)) < 0.5 and abs(Rnc - (HHc + LEc)) < 0.5: break

    return  Rn, G, HH, LE, HHg, LEg, HHc, LEc, Tg, Tc, count1


###############################################################################
# Subroutines 2003/12/10  by Atsushi Maruyama
#
###############################################################################
#--------------------------------------------------------------------------
#           Caliculation of Bulk transfer coefficients
#--------------------------------------------------------------------------
def cal_CX(z, d, zX1, zX2, psi1, psi2):

    if zX1 == 0 or zX2 == 0:
        return 0
    else:
        return Karman ** 2 * (log(z - d) - log(zX1) + psi1) ** -1 * (log(z - d) - log(zX2) + psi2) ** -1
        

#--------------------------------------------------------------------------
#           Calculation of aerodynamical parameters
#--------------------------------------------------------------------------
def cal_AeroPara(h, LAI, z0s, zTs, zqs, cdf, chf, cef):
    Aplus = cdf * LAI / float(2 * Karman ** 2)
    d = h * (1 - (1 - exp(-Aplus)) / Aplus)

    Term_z0 = (1 - exp(-Aplus) + exp(-2 * Aplus) * (-log(z0s / h)) ** (-1 / 0.45)) ** 0.45
    
    Term_z0t = cal_Term_zXt(z0s, h, Aplus)
    Term_zTt = cal_Term_zXt(zTs, h, Aplus)
    Term_zqt = cal_Term_zXt(zqs, h, Aplus)
    
    FT = chf / float(cdf)
    Fq = cef / float(cdf)
    
    Term_zT = cal_Term_zX(Term_z0, Term_zTt, FT, Aplus)
    Term_zq = cal_Term_zX(Term_z0, Term_zqt, Fq, Aplus)

    z0 = (h - d) * exp(-1 * Term_z0 ** -1)
    zT = (h - d) * exp(-1 * Term_zT ** -1)
    zq = (h - d) * exp(-1 * Term_zq ** -1)
    z0t = (h - d) * exp(-1 * Term_z0t ** -1)
    zTt = (h - d) * exp(-1 * Term_zTt ** -1)
    zqt = (h - d) * exp(-1 * Term_zqt ** -1)

    return d, z0, zT, zq, z0t, zTt, zqt


def cal_Term_zXt(zXs, h, Aplus):
    """Caliculate the values of z0t,zTt,zqt (z0,zT,zq when Fx=0)"""
    P0 = zXs / h
    P1 = 0.0115 * P0 ** 0.1 * exp(5 * P0 ** 0.22)
    P2 = 0.55 * exp(-0.58 * P0 ** 0.35)
    return -1 / (log(P0)) * (P1 / (P1 + Aplus * exp(Aplus))) ** P2

def cal_Term_zX(Term_z0, Term_zXt, FX, Aplus):
    """Caliculate the values of zT and zq"""
    CX0 = Term_z0 * Term_zXt
    CXm = (sqrt(1 + 8 * FX) - 1) / 2.0
    P3 = (FX + 0.084 * exp(-15 * FX)) ** 0.15
    P4 = 2 * FX ** 1.1
    TermxTerm = CXm * (1 - exp(-P3 * Aplus) + (CX0 / CXm) ** (1 / 0.9) * exp(-P4 * Aplus)) ** 0.9
    return TermxTerm / Term_z0

def cal_Ucmean(z, h, LAI, cdf, Uz):
    """Calculation of mean wind speed in the canopy(assuming exponential wind profile)"""
    # Dim gamma As Double
    # Dim Uh As Double        #Wind speed at the top of canopy(m/s)
        
    gamma = cdf * (LAI / h) / float(2 * Karman ** 2)
    Uh = Uz / (1 + log(gamma * (z - h) + 1))
    return Uh / (gamma * h) * (1 - exp(-gamma * h))

#--------------------------------------------------------------------------
#Caliculation of correction terms expressing the thermal stability effect
#--------------------------------------------------------------------------
def cal_psi(theta, CM, U, HH, cp, rho, z, d, z0, zT, zq):
    # Calculation of MoninL (Monin-Obukhov Length)
    ustar = sqrt(CM) * U
    if HH == 0: HH = 0.000001
    MoninL = -((theta + 273.15) * ustar ** 3) / (Karman * Gravitation * (HH / (cp * rho)))
    if MoninL == 0: MoninL = 0.000001
    
    # Calculation of zeta (Non-dimensional height)
    zeta = (z - d) / MoninL
    zeta0 = z0 / MoninL
    zetaT = zT / MoninL
    zetaq = zq / MoninL

    # Calculation of psi (correction terms that express the thermal stability effect)
    if zeta < 0:
        # Under unstable condtion
        x = (1 - 16 * zeta) ** (1 / 4.0)
        x0 = (1 - 16 * zeta0) ** (1 / 4.0)
        y = (1 - 16 * zeta) ** (1 / 2.0)
        y0 = (1 - 16 * zeta0) ** (1 / 2.0)
        yT = (1 - 16 * zetaT) ** (1 / 2.0)
        yq = (1 - 16 * zetaq) ** (1 / 2.0)
        
        psiM = log(((x0 ** 2 + 1) * (x0 + 1) ** 2) / ((x ** 2 + 1) * (x + 1) ** 2)) + 2 * (atan(x) - atan(x0))
        psiH = 2 * log((yT + 1) / (y + 1))
        psiE = 2 * log((yq + 1) / (y + 1))
        
    else:
        # Under stable condition
        alphaM = 3
        betaM = 10
        gammaM = 7 / 3.0
        alphaH = 7 / 400.0
        betaH = 0.005
        gammaH = 400

        psiM = gammaM * log((1 + alphaM * zeta + betaM * zeta ** 3) / (1 + alphaM * zeta0 + betaM * zeta0 ** 3))
        psiH = gammaH * log((1 + alphaH * zeta + betaH * zeta ** 2) / (1 + alphaH * zetaT + betaH * zetaT ** 2))
        psiE = gammaH * log((1 + alphaH * zeta + betaH * zeta ** 2) / (1 + alphaH * zetaq + betaH * zetaq ** 2))
    return  psiM, psiH, psiE


#########################################################################################
# Subroutine to calculate roughness length 2017/08/01 by Atsushi Maruyama
# (Maruyama et al. [2017] Water Resour. Res. 53, doi:10.1002/2017WR021019)

def cal_Roughness(z, h_geo, ustar_geo):
    KV = 1.53e-5            # Kinematic viscosity of the air at 20C (m2/s)
    TD = 2.15e-5            # Thermal diffusivity of the air at 20C (m2/s)

    Reynolds_h = ustar_geo * h_geo / KV
    
    if Reynolds_h < 16:
        Rfunc = (1 / Karman) * log(Reynolds_h) + 5.5
    else:
        Rfunc = 4
    
    InvTerm_z0s = (1 / Karman) * log(z / h_geo) + Rfunc
    z0s = z * exp(-1 * Karman * InvTerm_z0s)

    Reynolds_z0 = ustar_geo * z0s / KV
    Prandtl = KV / TD
    
    if Reynolds_z0 < 2:
        InvBH = -2.7
    else:
        InvBH = 8.5 * (Reynolds_z0 ** 0.25) * (Prandtl ** 0.5) - 4
    
    InvTerm_zTs = InvTerm_z0s + InvBH
    zTs = z * exp(-1 * Karman * InvTerm_zTs)
    zqs = zTs
    return z0s, zTs, zqs


#########################################################################################
# Subroutine of radiation parameters 2006/02/08 by Atsushi Maruyama
# (Maruyama et al. [2007] Jpn.Agric.Res.Quart.41,39-45)

def cal_RadPara(LAI, F, alpha, refg):
    tauS = exp(-sqrt(alpha) * F * LAI)                  # Canopy transmissivity for short-wave radiation
    tauL = exp(-F * LAI)                                # Canopy transmissivity for long-wave radiation
    refmax = ((1 - sqrt(alpha)) / (1 + sqrt(alpha)))    # Asymptotic value of albedo on the canopy
    ref = refmax - (refmax - refg) * (tauS ** 2)        # Albedo on the canopy
    return ref, tauS, tauL

def cal_F(DVS):
    """Calclation of F from developmental stage (DVS)
       F: Leaf inclination factor (extinction coefficient for black body)"""
    F0 = 0.21                 # F at transplanting
    F1 = 0.52                 # F at heading
    F2 = 1.1                  # F at maturity
    F_theta1 = 1.78           # Parameter expressing the increase after transplanting
    F_theta2 = 1.91           # Parameter expressing the increase after heading

    return F1 + (F0 - F1) * exp(-1 * (F_theta1 * DVS) ** 2) + (F2 - F1) * exp(-1 * (F_theta2 * (2 - DVS)) ** 2)

#########################################################################################
# Subroutine of bulk stomatal conductance 2006/03/21 by Atsushi Maruyama
# (Maruyama and Kuwagata [2008] Agric.For.Met.148,1161-1173)

def cal_gsf(DVS, Scmean):
    """Calclation of gsf from developmental stage (DVS) and scmean
       gsf: bulk stomatal conductance (m/s)
       Scmean: mean short-wave radiation absorbed by the canopy (W/m2)"""
    # Dim gsmax As Double     # maximum bulk stomatal conductance (m/s)
    gsmax0 = 0.06             # gsmax at transplanting (m/s)
    gsmax1 = 0.023            # gsmax at heading (m/s)
    gsmax2 = 0.009            # gsmax at maturity (m/s)
    gs_theta1 = 1.65          # Parameter expressing the decrease after transplanting
    gs_theta2 = 1.08          # Parameter expressing the decrease after heading
    Scmean05 = 66             # Scmean when gs value is half of gsmax (W/m2)
    
    gsmax = gsmax1 + (gsmax0 - gsmax1) * exp(-1 * (gs_theta1 * DVS) ** 8) + (gsmax2 - gsmax1) * exp(-1 * (gs_theta2 * (2 - DVS)) ** 8)
    return gsmax * Scmean / (Scmean05 + Scmean)


##############################################################################################
# Sigmoid function to calculate LAI from DVS
# (Maruyama [2017] not published)

def cal_LAI(DVS, LAImax):

    DVS_h1 = 0.65       # DVS when LAI value is half of LAImax after transplanting
    DVS_h2 = 1.84       # DVS when LAI value is half of LAImax after heading
    LAI_K1 = 7.6        # Parameter affecting the degree of LAI increase after transplanting
    LAI_K2 = 2.4        # Parameter affecting the degree of LAI decrease after heading

    return LAImax * (1 - 0.5 ** exp(-1 * LAI_K1 * (DVS_h1 - DVS)) - 0.5 ** exp(-1 * LAI_K2 * (DVS - DVS_h2)))


##############################################################################################
# Allometric model to calculate canopy height from LAI
# (Maruyama and Kuwagata [2010] Agric.For.Met.150,919-930)

def cal_h(DVS, LAI):

    Aro1 = 0.458
    Bro1 = 0.367

    Aro2 = 0.714
    Bro2 = 0.285

    if DVS < 1:
        return Aro1 * LAI ** Bro1
    else:
        return Aro2 * LAI ** Bro2

##############################################################################################

def DSMcalculation(TG, LAT,ALT,DOY,WIND,TMP_mea,TMP_max,TMP_min,Sday,_DVS=0,LAImax=4.6,Dw=50,refg=0.1):

    eave = cal_e(TMP_min, 100)

    (Tgpast, Rn_ave,G_ave,HH_ave,LE_ave,Tg_ave,Tc_ave,Tg_max,Tc_max,Tg_min,Tc_min
    ) = DSMcalculation_impl(LAT,ALT,DOY,WIND,TMP_mea,TMP_max,TMP_min,eave,Sday,_DVS,LAImax,Dw,refg, TG)

    return Tgpast, Tg_ave,Tg_max,Tg_min

# same as VBA I/O
def DSMcalculation_impl(LAT,ALT,DOY,WIND,TMP_mea,TMP_max,TMP_min,eave,Sday,_DVS,LAImax,Dw,refg, TG):

    #LAImax = 4.6  # Maximum laef area index (m2/m2)

    TgM = 0
    Tg0 = 0
    
    if None in TG:
        Tgpast = 25 * [0] 
        Spinup = 5
    else:
        Tgpast = TG
        Spinup = 0

    za = 2

    Latitude = LAT * pi / 180
    Elevation = ALT
    Uave = WIND
    Tave = TMP_mea
    Tmax = TMP_max
    Tmin = TMP_min
    _LAI = cal_LAI(_DVS, LAImax)
    _h = cal_h(_DVS, _LAI)

    Patm = cal_Patm(Elevation)
    N0, Lday = cal_Lday_fromS(DOY, Latitude, Tave, eave, Patm, Sday)

    TgM = sum(Tgpast[1:])/24 
    Tg0 = Tgpast[0]

    while (Spinup >= 0):
        #-------Reset of daily mean values---------
        Rn_ave = 0
        G_ave = 0
        HH_ave = 0
        LE_ave = 0
        HHg_ave = 0
        LEg_ave = 0
        HHc_ave = 0
        LEc_ave = 0
        Tg_ave = 0
        Tc_ave = 0
        Tg_max = -99
        Tc_max = -99
        Tg_min = 99
        Tc_min = 99
        #-----------------------------------------
        
        for TimeCount in range(1,25):
  
            Ta = cal_Ta(Tave, Tmax, Tmin, TimeCount)
            ea = eave
            Ua = cal_Ua(Uave, TimeCount)    # 2012/12 add diurnal variation in U
            Sd = cal_Sd(Sday, TimeCount)
            Ld = Lday

            (Rn, G, HH, LE, HHg, LEg, HHc, LEc, Tg, Tc, count1
            ) = DSM(Patm, Ta, ea, Ua, Sd, Ld, _h, _LAI, _DVS, za, Tg0, TgM, Dw, refg)

            #----Calculation of new TgM and Tg0 values----
            Tgpast = [Tg] + Tgpast[:24] 
            TgM = sum(Tgpast[1:])/24 
            Tg0 = Tg

            #-------Calculation of daily mean values-------
            Rn_ave += Rn / 24
            G_ave += G / 24
            HH_ave += HH / 24
            LE_ave += LE / 24
            HHg_ave += HHg / 24
            LEg_ave += LEg / 24
            HHc_ave += HHc / 24
            LEc_ave += LEc / 24
            Tg_ave += Tg / 24
            Tc_ave += Tc / 24
            Tg_max = max(Tg_max, Tg)
            Tc_max = max(Tc_max, Tc)
            Tg_min = min(Tg_min, Tg)
            Tc_min = min(Tc_min, Tc)
            #-----------------------------------------------

        Spinup -= 1
    return Tgpast, Rn_ave,G_ave,HH_ave,LE_ave,Tg_ave,Tc_ave,Tg_max,Tc_max,Tg_min,Tc_min


def WaterTemp(timedomain, lalodomain, dvs, Dw, laim=4.6, alb=0.1, ws=None, cli=False):
    assert ((cli == False) or (ws != None)), "Error : please set ws for cli=True"
    ALTs  = AMD.GetGeoData("altitude", lalodomain)
    GSRs  = AMD.GetMetData("GSR",  timedomain, lalodomain, cli=cli)
    tim = GSRs[1]
    lat = GSRs[2]
    lon = GSRs[3]
    assert (len(tim) == len(dvs)), "Error : date list and dvs list have different length"
    assert (len(tim) == len(Dw)), "Error : date list and Dw list have different length"
    TMP_maxs = AMD.GetMetData("TMP_max", timedomain, lalodomain, cli=cli)
    TMP_mins = AMD.GetMetData("TMP_min", timedomain, lalodomain, cli=cli)
    TMP_meas = AMD.GetMetData("TMP_mea", timedomain, lalodomain, cli=cli)

    if ws is None:
        WINDs = AMD.GetMetData("WIND", timedomain, lalodomain, cli=cli)
    else:
        WIND = ws

    Tmea = np.zeros((len(tim), len(lat), len(lon)), dtype=np.float32); Tmea[:,:,:] = np.nan
    Tmax = np.zeros((len(tim), len(lat), len(lon)), dtype=np.float32); Tmax[:,:,:] = np.nan
    Tmin = np.zeros((len(tim), len(lat), len(lon)), dtype=np.float32); Tmin[:,:,:] = np.nan
    for la in range(len(lat)):
        print(la+1,"/",len(lat))	#進捗が分かったらいいかなと思って。-20170131大野追加
        for lo in range(len(lon)):
            TG = 25*[None]
            ALT = ALTs[0][la,lo]
            if np.isnan(ALT): continue
            for i in range(len(tim)):
                skipCalc = False
                if ws is None:
                    WIND = WINDs[0][i,la,lo]
                if np.isnan(WIND):
                    skipCalc = True
                TMP_mea = TMP_meas[0][i,la,lo]
                if np.isnan(TMP_mea):
                    skipCalc = True
                TMP_max = TMP_maxs[0][i,la,lo]
                if np.isnan(TMP_max):
                    skipCalc = True
                TMP_min = TMP_mins[0][i,la,lo]
                if np.isnan(TMP_min):
                    skipCalc = True
                if np.isnan(GSRs[0][i,la,lo]):
                    skipCalc = True
                else:
                    GSR = GSRs[0][i,la,lo]*GSR_FACTOR
                DOY = (tim[i]-dt(tim[i].year,1,1)).days + 1
                if skipCalc:
                    TG = 25*[None]
                else:
                    ret = DSMcalculation(TG, lat[la], ALT, DOY, WIND,
                                         TMP_mea, TMP_max, TMP_min,
                                         GSR, dvs[i], laim, Dw[i], alb)
                    TG = ret[0]
                    Tmea[i,la,lo] = ret[1]
                    Tmax[i,la,lo] = ret[2]
                    Tmin[i,la,lo] = ret[3]
    return Tmea.astype(np.float32), Tmax.astype(np.float32), Tmin.astype(np.float32), \
           tim, lat, lon
