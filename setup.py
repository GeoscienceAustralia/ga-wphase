from numpy.distutils.core import Extension, setup
import os

def get_package_data(base_dir, sub_dir, base_items):
    cdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), base_dir, sub_dir)
    skip = len(cdir) - len(sub_dir)

    for dp, dn, fns in os.walk(cdir):
        dp = dp[skip:]
        for fn in fns:
            base_items.append('{}/{}'.format(dp, fn))

    return base_items

wphase_files = get_package_data(
    'wphase', 'ATWS', [])

# This closed-source fortran code may or may not be present.
# Only try to compile it if it's there.
FORTRAN_BPFILTER = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
        'wphase', 'psi', 'bpfilter.f')

if os.path.exists(FORTRAN_BPFILTER):
    extra_extensions = [
        Extension('wphase.psi.bpfilter', ['wphase/psi/bpfilter.f']),
    ]
else:
    extra_extensions = []

setup(
    name = 'Wphase',
    version = '0.1',
    packages = ['wphase', 'wphase.psi'],
    author = 'Geoscience Austraila',
    description = 'Wphase calculations and web interface.',
    package_dir = {'wphase':'wphase', 'wphase.psi': 'wphase/psi'},
    scripts = ['scripts/wphase'],
    package_data = {
        'wphase': wphase_files,
        'wphase.psi': ['tables/*.*']},
    ext_modules = extra_extensions + [
        Extension('wphase.psi._libtau', [
            'wphase/psi/_libtau.pyf',
            'wphase/psi/ttimes_subrout.f',
            'wphase/psi/libtau.f'])
    ],
    install_requires = ['numpy', 'pytest', 'obspy', 'matplotlib', 'scipy>=0.16.0', 'pyinstrument==0.13.1']
    )
