import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

def make_spot(nside, spot_theta, spot_phi, spot_radius, spot_intensity):
    """
    スポットが1つあるマップを作る。
    
    Parameters
    ----------
    nside : int
        マップの解像度を決めるパラメータで、2の累乗。
    spot_theta : float
        スポット中心のθ座標(rad)。
    spot_phi : float
        スポット中心のφ座標(rad)。
    spot_radius : float
        スポットの半径(rad)。
    spot_intensity : float
        スポットの強度(0が真っ黒)。
    
    Examples
    ----------
    map = make_spot(np.pi/2, np.pi/4, np.pi/10, 0.2)
    hp.mollview(map, title="Star Surface with Spot", cmap="inferno", unit="Intensity")
    plt.savefig('spot.png')
    """

    # 星表面の初期条件
    npix = hp.nside2npix(nside)  # ピクセル数の計算
    map = np.ones(npix)  # 均一な輝度を設定

    # 各ピクセルの座標（θ, φ）を取得
    theta, phi = hp.pix2ang(nside, np.arange(npix))

    # スポットの中心から各ピクセルまでの角距離を計算
    angular_distance = hp.rotator.angdist([spot_theta, spot_phi], [theta, phi])

    # スポットの範囲内にあるピクセルを見つける
    spot_mask = angular_distance < spot_radius

    # スポットを強度に応じて設定
    map[spot_mask] = spot_intensity

    return map

def add_spot(map, spot_theta, spot_phi, spot_radius, spot_intensity):
    """
    マップにスポットを1つ追加する。
    
    Parameters
    ----------
    map : numpy.ndarray
        HEALPixによる輝度マップ(長さがピクセル数の1次元配列)。
    spot_theta : float
        スポット中心のθ座標(rad)。
    spot_phi : float
        スポット中心のφ座標(rad)。
    spot_radius : float
        スポットの半径(rad)。
    spot_intensity : float
        スポットの強度(0が真っ黒)。
    
    Examples
    ----------
    map = np.ones(hp.nside2npix(128))
    map = add_spot(map, np.pi/3, np.pi/4, np.pi/10, 0.2) #左目
    map = add_spot(map, np.pi/2.5, np.pi/4, np.pi/10, 0.2) #左目
    map = add_spot(map, np.pi/3, -np.pi/4, np.pi/10, 0.2) #右目
    map = add_spot(map, np.pi/2.5, -np.pi/4, np.pi/10, 0.2) #右目
    map = add_spot(map, np.pi/1.4, 0, np.pi/10, 0.5) #口
    map = add_spot(map, np.pi/1.7, np.pi/2, np.pi/10, 0.8) #頬
    map = add_spot(map, np.pi/1.7, -np.pi/2, np.pi/10, 0.8) #頬
    hp.mollview(map, title="Star Surface with Spot", min=0., cmap="inferno", unit="Intensity")
    plt.savefig('spot.png')
    """

    # 星表面の初期条件
    nside = int(np.sqrt(len(map)/12))  # Nsideは2の累乗
    npix = len(map)

    # 各ピクセルの座標（θ, φ）を取得
    theta, phi = hp.pix2ang(nside, np.arange(npix))

    # スポットの中心から各ピクセルまでの角距離を計算
    angular_distance = hp.rotator.angdist([spot_theta, spot_phi], [theta, phi])

    # スポットの範囲内にあるピクセルを見つける
    spot_mask = angular_distance < spot_radius

    # スポットを強度に応じて設定
    map[spot_mask] *= spot_intensity

    return map
