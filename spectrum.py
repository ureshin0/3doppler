import healpy as hp
import numpy as np

from const import Const

def incline(theta0, phi0, alpha):
    # Convert spherical to Cartesian
    x0 = np.sin(theta0) * np.cos(phi0)
    y0 = np.sin(theta0) * np.sin(phi0)
    z0 = np.cos(theta0)

    # Apply rotation matrix around the y-axis
    x = np.cos(alpha) * x0 + np.sin(alpha) * z0
    y = y0
    z = -np.sin(alpha) * x0 + np.cos(alpha) * z0

    # Convert back to spherical
    theta = np.arccos(z)
    phi = np.arctan2(y, x)

    return theta, phi

def limb_darkening(u, mu):
    return 1-u*(1-mu)

def doppler_shift(vlos):
    """
    視線速度からドップラーシフトの倍率を計算する。

    Parameters
    ----------
    v : float
        視線速度(km/s)。視線から遠ざかる向きが正。
    """
    c = Const().c0 * 10**(-3)  # 光速 (km/s)
    beta = vlos/c
    return (1 + beta) / np.sqrt(1 - beta**2)

def observe_spectrum(map, vrot, wavelengths, line_profile, inclination=np.pi/2, phase=0., u=0, output_resolution=1000 ,normalize=True):
    """
    マップを元に、ドップラー効果を考慮したスペクトルを生成する。

    Parameters
    ----------
    map : array
        星表面のマップ。
    vrot : float
        星の自転速度(km/s)。
    wavelengths : array
        一様なスペクトルのシフト前の波長(nm)。
    line_profile : array
        一様なスペクトル強度。
    inclination : float
        自転軸の傾き(rad)。pi/2のとき視線に垂直。
    phase : float
        観測するフェーズ。1.0で1自転。
    u : float
        周縁減光の強さを表す係数。0から1まで。
    output_resolution : int
        出力スペクトルの解像度。
    normal : boolean
        最後に正規化するかどうか。
    """

    nside = hp.npix2nside(len(map))
    theta0, phi0 = hp.pix2ang(nside, np.arange(len(map)))

    # 自転
    phi0 += phase*2*np.pi

    # 観測者からの視線速度を計算
    vlos = vrot * np.cos(np.pi/2-inclination) * np.sin(theta0) * np.sin(phi0)

    # 観測者から見た座標に変換
    theta, phi = incline(theta0, phi0, np.pi/2-inclination)

    # スペクトル波長範囲を設定
    wl_min = wavelengths[0] * doppler_shift(-vrot)
    wl_max = wavelengths[-1] * doppler_shift(vrot)
    observed_wavelengths = np.linspace(wl_min, wl_max, output_resolution)

    # 観測されるスペクトルを初期化
    observed_spectrum = np.zeros_like(wavelengths)

    # 各ピクセルについて放射をドップラーシフトさせ、スペクトルに足し合わせる
    for i in range(len(map)):
        if -np.pi/2 < phi[i] < np.pi/2:
            pixel_intensity = map[i]
            shifted_wavelengths = wavelengths * doppler_shift(vlos[i])
            mu = np.sin(theta[i]) * np.cos(phi[i])
            # シフトされた波長に基づいてスペクトルに追加
            observed_spectrum += pixel_intensity * np.interp(observed_wavelengths, shifted_wavelengths, line_profile) * limb_darkening(u,mu) * np.sin(theta[i]) * np.cos(phi[i])

    # スペクトルの正規化
    if normalize:
        observed_spectrum /= np.max(observed_spectrum)

    return observed_wavelengths, observed_spectrum


if __name__ == "__main__":
    from map import make_spot
    from scipy.stats import norm
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt

    # 生成例
    spot_theta = np.pi / 2
    spot_phi = -np.pi / 2
    spot_radius = np.pi / 4
    spot_intensity = 0.1

    vrot = 20  # 回転速度 (km/s)
    wl0=656.28
    wl = np.linspace(wl0-0.04, wl0+0.04, 1000)
    line_profile = norm.pdf(wl, wl0, 0.01)  # 波長0.01nm幅のガウシアン

    map = make_spot(64, spot_theta, spot_phi, spot_radius, spot_intensity)

    hp.mollview(map, title="Star Surface with Spot", min=0., cmap="inferno", unit="Intensity", flip='geo')
    plt.savefig('spot.png')
    plt.clf()

    for i in range(8):
        wavelengths, spectrum = observe_spectrum(map, vrot, wl, line_profile, phase=0.125*i, normalize=False)
        plt.plot(wavelengths, spectrum, label=f"t={i}T/8", color=cm.hot(i/8), zorder=8-i)
    plt.legend()
    plt.title("i=90deg")
    plt.savefig('spectrum_rotation_90.png')
    plt.clf()

    for i in range(8):
        wavelengths, spectrum = observe_spectrum(map, vrot, wl, line_profile, inclination=np.pi/4, phase=0.125*i, normalize=False)
        plt.plot(wavelengths, spectrum, label=f"t={i}T/8", color=cm.hot(i/8), zorder=8-i)
    plt.legend()
    plt.title("i=45deg")
    plt.savefig('spectrum_rotation_45.png')
    plt.clf()

    for i in range(8):
        wavelengths, spectrum = observe_spectrum(map, vrot, wl, line_profile, inclination=0, phase=0.125*i, normalize=False)
        plt.plot(wavelengths, spectrum, label=f"t={i}T/8", color=cm.hot(i/8), zorder=8-i)
    plt.legend()
    plt.title("i=0deg")
    plt.savefig('spectrum_rotation_0.png')
    plt.clf()
