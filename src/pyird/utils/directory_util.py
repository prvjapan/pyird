import shutil


def cp(anadir, fitsset, iraftag='ap'):
    fitsset_path = fitsset.path(string=False, check=True)[0]
    fitsset_path_ap = fitsset.fitsdir/'database' / \
        (iraftag+fitsset.tag+fitsset.extension)
    current_fitsset_path = anadir/(fitsset.tag+fitsset.extension+'.fits')
    current_fitsset_path_ap = anadir/'database' / \
        (iraftag+fitsset.tag+fitsset.extension)
    if not current_fitsset_path.exists():
        print('COPYING:', fitsset_path, current_fitsset_path)
        shutil.copyfile(fitsset_path, current_fitsset_path)

    if not current_fitsset_path_ap.exists():
        print('COPYING:', fitsset_path_ap, current_fitsset_path_ap)

        shutil.copyfile(fitsset_path_ap, current_fitsset_path_ap)
    return current_fitsset_path
