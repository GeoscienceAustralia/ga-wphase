# -*- coding: utf-8 -*-
"""
Modified version of the old Obspy taup implementation.
The new version is very slow for near real time apllications.
In order to install it you must go to wphase.psi and run:

f2py -c _libtau.pyf ttimes_subrout.f libtau.f
"""

import os
import wphase.psi._libtau as _libtau

ttimes = _libtau.ttimes
_taup_dir = os.path.dirname(os.path.abspath(__file__))




AVAILABLE_PHASES = [
    'P', "P'P'ab", "P'P'bc", "P'P'df", 'PKKPab', 'PKKPbc', 'PKKPdf', 'PKKSab',
    'PKKSbc', 'PKKSdf', 'PKPab', 'PKPbc', 'PKPdf', 'PKPdiff', 'PKSab', 'PKSbc',
    'PKSdf', 'PKiKP', 'PP', 'PS', 'PcP', 'PcS', 'Pdiff', 'Pn', 'PnPn', 'PnS',
    'S', "S'S'ac", "S'S'df", 'SKKPab', 'SKKPbc', 'SKKPdf', 'SKKSac', 'SKKSdf',
    'SKPab', 'SKPbc', 'SKPdf', 'SKSac', 'SKSdf', 'SKiKP', 'SP', 'SPg', 'SPn',
    'SS', 'ScP', 'ScS', 'Sdiff', 'Sn', 'SnSn', 'pP', 'pPKPab', 'pPKPbc',
    'pPKPdf', 'pPKPdiff', 'pPKiKP', 'pPdiff', 'pPn', 'pS', 'pSKSac', 'pSKSdf',
    'pSdiff', 'sP', 'sPKPab', 'sPKPbc', 'sPKPdf', 'sPKPdiff', 'sPKiKP', 'sPb',
    'sPdiff', 'sPg', 'sPn', 'sS', 'sSKSac', 'sSKSdf', 'sSdiff', 'sSn']


def getTravelTimes(delta, depth, model='iasp91'):
    """
    Returns the travel times calculated by the iaspei-tau, a travel time
    library by Arthur Snoke (http://www.iris.edu/pub/programs/iaspei-tau/).

    :type delta: float
    :param delta: Distance in degrees.
    :type depth: float
    :param depth: Depth in kilometer.
    :type model: str, optional
    :param model: Either ``'iasp91'`` or ``'ak135'`` velocity model. Defaults
        to ``'iasp91'``.
    :rtype: list of dicts
    :return:
        A list of phase arrivals given in time order. Each phase is represented
        by a dictionary containing phase name, travel time in seconds, take-off
        angle, and various derivatives (travel time with respect to distance,
        travel time with respect to depth and the second derivative of travel
        time with respect to distance).

    .. rubric:: Example

    >>> from obspy.taup.taup import getTravelTimes
    >>> tt = getTravelTimes(delta=52.474, depth=611.0, model='ak135')
    >>> len(tt)
    24
    >>> tt[0]  #doctest: +SKIP
    {'phase_name': 'P', 'dT/dD': 7.1050525, 'take-off angle': 45.169445,
     'time': 497.53741, 'd2T/dD2': -0.0044748308, 'dT/dh': -0.070258446}
    """
    # Raise an error, otherwise libtau sends an EXIT signal. Depends on the
    # model but 800 km works for the included models.
    if depth > 800.00:
        raise ValueError("Source depth of %.2f km is too deep." % depth)
    model_path = os.path.join(_taup_dir, 'tables', model)
    if not os.path.exists(model_path + os.path.extsep + 'hed') or \
       not os.path.exists(model_path + os.path.extsep + 'tbl'):
        msg = 'Model %s not found' % model
        raise ValueError(msg)

    # Depth in kilometer.
    depth = abs(depth)

    # modnam is a string with 500 chars.
    modnam = os.path.join(_taup_dir, 'tables', model).ljust(500)

    phase_names, tt, toang, dtdd, dtdh, dddp = ttimes(delta, depth, modnam)

    phases = []
    for _i, phase in enumerate(phase_names):
        # An empty returned string will contain "\x00".
        phase_name = phase.tostring().strip().\
            replace(b"\x00", b"").decode()
        if not phase_name:
            break
        time_dict = {
            'phase_name': phase_name,
            'time': tt[_i],
            'take-off angle': toang[_i],
            'dT/dD': dtdd[_i],
            'dT/dh': dtdh[_i],
            'd2T/dD2': dddp[_i]}
        phases.append(time_dict)
    return phases

def getPtime(delta, depth, model='iasp91'):
    """
    Returns the travel times calculated by the iaspei-tau, a travel time
    library by Arthur Snoke (http://www.iris.edu/pub/programs/iaspei-tau/).

    :type delta: float
    :param delta: Distance in degrees.
    :type depth: float
    :param depth: Depth in kilometer.
    :type model: str, optional
    :param model: Either ``'iasp91'`` or ``'ak135'`` velocity model. Defaults
        to ``'iasp91'``.
    :rtype: list of dicts
    :return: P arrival time in seconds.



    .. rubric:: Example

    >>> from obspy.taup.taup import getTravelTimes
    >>> t_p = getPtime(delta=52.474, depth=611.0, model='ak135')
    497.537
    """
    # Raise an error, otherwise libtau sends an EXIT signal. Depends on the
    # model but 800 km works for the included models.
    if depth > 800.00:
        raise ValueError("Source depth of %.2f km is too deep." % depth)
    model_path = os.path.join(_taup_dir, 'tables', model)
    if not os.path.exists(model_path + os.path.extsep + 'hed') or \
       not os.path.exists(model_path + os.path.extsep + 'tbl'):
        msg = 'Model %s not found' % model
        raise ValueError(msg)

    # Depth in kilometer.
    depth = abs(depth)

    # modnam is a string with 500 chars.
    modnam = os.path.join(_taup_dir, 'tables', model).ljust(500)

    phase_names, tt, toang, dtdd, dtdh, dddp = ttimes(delta, depth, modnam)

    return tt[0]




if __name__ == '__main__':
    import doctest
    doctest.testmod(exclude_empty=True)
