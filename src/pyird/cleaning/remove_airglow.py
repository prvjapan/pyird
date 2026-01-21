import pandas as pd
from importlib import resources
import matplotlib.pyplot as plt

def load_oliva15():
    """
    Returns:
        all wavlengths of airglow lines in nm
    """
    file_tb1 = resources.files('pyird').joinpath('data/Oliva+15_table1.dat')
    file_tb2 = resources.files('pyird').joinpath('data/Oliva+15_table2.dat')
    oliva_tb1 = pd.read_csv(file_tb1, sep='\s+', names=['lambda1', 'Iden1', 'lambda2', 'Iden2', 'Flux'])
    oliva_tb2 = pd.read_csv(file_tb2, sep='\s+', names=['lambda', 'Iden', 'Flux', 'q_Flux'])
    wav_oh = pd.concat([oliva_tb1['lambda1'], oliva_tb1['lambda2']])
    wav_all = pd.concat([oliva_tb2['lambda'], wav_oh])
    wav_all = wav_all.sort_values(ignore_index=True)
    wav_all = wav_all * 1e-1 #AA to nm
    return wav_all


def wav_around_airglow(wav_obs, wav_airglow, mask_width=0.035):
    """
    Args:
        wav_obs: wavelengths for observed data
        wav_airglow: wavlengths of airglow mask
        mask_width: masking width in nm
    
    Returns:
        index for wav_obs that matches {wavelengths in wav_airglow} +/- {mask_width}

    Notes:
        O2 bandhead features cannot masked by this method.
    """
    ind_mask = []
    for wav_airglow_i in wav_airglow:
        ind_mask_i = wav_obs[(wav_airglow_i-mask_width <= wav_obs) & (wav_obs <= wav_airglow_i+mask_width)].index.values
        ind_mask.extend(ind_mask_i)
    return ind_mask

def df_mask_airglow(df_obs, ind_mask, plot=False):
    """
    Args:
        df_obs: pandas DataFrame 
        ind_mask: index for masking
        plot: if True, plot before and after masked spectra

    Returns:
        masked DataFrame
    """
    df_masked = df_obs.drop(ind_mask)
    if plot:
        fig, ax = plt.subplots()
        ax.plot(df_obs["wav"], df_obs["flux"], "k", lw=1, label="obs")
        ax.plot(df_masked["wav"], df_masked["flux"], "tab:blue", lw=1, label="masked")
        ax.legend()
        #ax.set(xlim=(1500, 1580), ylim=(0,0.02))
        plt.show()
    return df_masked


if __name__ == '__main__':
    date_dir = "202105"
    date_file = "202105"
    target = "GD165B"
    band = "h"
    save = False

    file = f"/Users/yuikasagi/IRD/PhDwork/pyird/data/{date_dir}/reduc_v1.1/ncw{target}_{date_file}_{band}_m2.dat"
    df_obs = pd.read_csv(file, sep="\s+", names=["wav","flux","uncertainty"])

    wav_airglow = load_oliva15()
    ind_mask = wav_around_airglow(df_obs["wav"], wav_airglow)
    df_masked = df_mask_airglow(df_obs, ind_mask, plot=True)

    save_path = f"/Users/yuikasagi/IRD/PhDwork/pyird/data/{date_dir}/reduc_v1.1/ncw{target}_{date_file}_{band}_m2_ohmasked.dat"
    if save:
        df_masked.to_csv(save_path, sep=" ", index=False, header=False)

