import matplotlib.pyplot as plt

def show_wavcal_spectrum(df_wspec,**kwargs):
    """show the wavelength calibrated spectrum.

    Args:
        df_wspec: pandas.DataFrame that contains wavelengths and flux
        kwargs: keys for plot

    """
    wav = df_wspec['wav'].values
    flux_col = df_wspec.columns[[('flux' in x) for x in df_wspec.columns]]
    flux = df_wspec[flux_col].values

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(wav,flux,**kwargs)
    plt.show()
