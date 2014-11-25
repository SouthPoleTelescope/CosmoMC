import chainDumpPlot as r
import CMBlikes
import os
import planckStyleNG
from pylab import *

def_figsize = 12
def_max = 400
diff_base = r'C:\tmp\Planck\final_Nov14\base_plikHM_TT_lowTEB.minimum.theory_cl'
commander = r'C:\Users\antony\Dropbox\Planck\final_Nov14\residuals\cls_commander_2014_rc_TT.dat'
plikLite = False
if plikLite:
    CMB_data_file = r'C:\Users\antony\Dropbox\Planck\final_Nov14\residuals\PlikLite_v17spectra.txt'
else:
    CMB_data_file = r'C:\Users\antony\Dropbox\Planck\final_Nov14\residuals\plik_v17_HM_Dl_CMBcoadded.dat'

lensing_dataset = r'C:\Work\Dist\git\cosmomcplanck\data\planck_lensing\smica_g30_ftl_full_pp_bandpowers.dat'

linePlot = True
make_file_dir = 'z://'
if make_file_dir:
    def_max = 100
    linePlot = False
    def_figsize = 3.5

def makeRainbow(ls, pars, cls, parname, lpow, delta_to, Lmax, color_pow=0.5, last_colors=None):
    if linePlot:
        return r.rainbowLinePlot(ls, pars, cls, parname=parname, lpow=lpow, delta_to=delta_to, alpha=1)
    else:
        return r.rainbowDensityPlot(ls, pars, cls, parname=parname, lpow=lpow, delta_to=delta_to, Lmax=Lmax, color_pow=color_pow, last_colors=last_colors)

def exportFile(chainRoot, parname):
    tight_layout()
    if make_file_dir:
        savefig(os.path.join(make_file_dir, os.path.split(chainRoot)[-1] + '_' + parname + '.pdf'), bbox_inches='tight', dpi=600)


def lensingRainbow(chainRoot, parname='omegach2', last_colors=None, Lmax=400, rescale=1e7):
    figure(figsize=(def_figsize, def_figsize * 3 / 4))
    cl_ix = 5
    ls, pars, cls = r.loadFiles(chainRoot, cl_ix=cl_ix, max_files=def_max, rescale=rescale)
    res = makeRainbow(ls, pars, cls, parname=parname, lpow=2, delta_to=None, Lmax=Lmax, last_colors=last_colors, color_pow=0.5)

    Lbin, pts, sig = CMBlikes.readTextCommentColumns(lensing_dataset, ['L_av', 'PP', 'Error'])
    pts *= rescale
    sig *= rescale
    errorbar(Lbin, pts, sig, fmt='.', color='k')

    ylabel(r'$[L(L+1)]^2C_L^{\phi\phi}/2\pi$\,[$\times 10^7$]')
#    setp(getp(cb.ax, 'ymajorticklabels'), fontsize=self.settings.colorbar_axes_fontsize)
    ylim([0, 1.8e-7 * rescale])
    xlim([2, Lmax])
    exportFile(chainRoot, parname)
    return res

def plotCommander(lmax=29):
    planck = np.loadtxt(commander)
    col = 'k'
    ls = planck[:lmax - 1, 0]
    print ls
    plot(ls, planck[:lmax - 1, 3], 'o' , color=col, ls='None', markersize=1)
    dat = planck
    for i in range(min(lmax - 1, planck.shape[0])):
        gca().add_line(Line2D((dat[i, 0], dat[i, 0]), (dat[i, 3] - dat[i, 5], dat[i, 3] + dat[i, 4]), color=col, linewidth=0.5))

def CMBRainbow(chainRoot, parname='omegach2', cl_ix=1, diff_to_file=None, Lmax=2000, data=True, ymax=None, tag=None):
    figure(figsize=(def_figsize, def_figsize * 3 / 4))
    calParam = 'calPlanck'
    if cl_ix == 2:
        lpow = 2
    else:
        lpow = 2
    ls, pars, cls = r.loadFiles(chainRoot, cl_ix=cl_ix, max_files=def_max, calParam=calParam)

    if diff_to_file: delta_to = loadtxt(diff_to_file)[:, cl_ix]
    else: delta_to = None
    res = makeRainbow(ls, pars, cls, parname=parname, lpow=lpow, delta_to=delta_to, Lmax=Lmax, color_pow=0.5)

    datapoints = loadtxt(CMB_data_file)
    Lbin = datapoints[:, 0]
    ioff = cl_ix * 2 - 1
    if not plikLite:
        if cl_ix == 2: ioff = 5
        if cl_ix == 3: ioff = 3
    pts = datapoints[:, ioff]
    sig = datapoints[:, ioff + 1]
    if lpow == 1: ylabel(r'$[\ell(\ell+1)]^{1/2} C^{TE}_\ell/2\pi$\,[$\mu{\rm K}^2$]')
    if lpow == 2 and cl_ix == 3: ylabel(r'$D^{EE}_\ell$\,[$\mu{\rm K}^2$]')
    if lpow == 2 and cl_ix == 1: ylabel(r'$D^{TT}_\ell$\,[$\mu{\rm K}^2$]')

    if diff_to_file:
        for i, L in enumerate(Lbin):
            ix = np.where(ls == L)
            pts[i] -= delta_to[ix]
    if data:
        sc = (Lbin * (Lbin + 1.)) ** (lpow / 2. - 1)
#        errorbar(Lbin, pts * sc, sig * sc, fmt='.', color='k', zorder=2, markersize=1, elinewidth=0.5, mew=0.5)
        col = 'k'
        plot(Lbin, pts * sc, 'o' , color=col, ls='None', markersize=1)
        dat = pts * sc
        for i, d in enumerate(dat):
            gca().add_line(Line2D((Lbin[i], Lbin[i]), (d - sc[i] * sig[i], d + sc[i] * sig[i]), color=col, linewidth=0.5))
        if cl_ix == 1:
            plotCommander()


    xlim([1, Lmax])
    xlabel(r'$\ell$')
#    if cl_ix == 2: ylim([-0.4, 0.4])
    if ymax: ylim([0, ymax])
    exportFile(chainRoot, tag or parname)
#    savefig(r'z:' + os.sep + os.path.split(chainRoot)[-1] + '_' + parname + '.pdf', bbox_inches='tight')
    return res

def lensingPaper():

    lensingRainbow(r'C:\tmp\Planck\final_Nov14\cl_chains\lensonly\base_lensonly_post', parname='omegach2')
    lensingRainbow(r'C:\tmp\Planck\final_Nov14\cl_chains\base_plikHM_TT_lowTEB_post', parname='omegach2')

if False:
    d = loadtxt(r'z:\base_plikHM_TT_lowTEB_lensing_lensedCls.dat')
    diso = loadtxt(r'z:\base_plikHM_TT_lowTEB_lensing_iso1_lensedCls.dat')
    dx = loadtxt(r'z:\base_plikHM_TT_lowTEB_lensing_isoxp1_lensedCls.dat')

    # print d.shape, diso.shape
    d = d[:2000, :]
    diso = diso[:2000, :]
    dx = dx[:2000, :]

    ls = d[:, 0]
    ix = 4
    subplot(211)
    plot(ls, d[:, ix], color='k')
    plot(ls, diso[:, ix], color='g')
    plot(ls, dx[:, ix], color='r')
    xlim([0, 1000])
    # gca().set_yscale('log')
    subplot(212)
    plot(ls, (dx[:, ix] - d[:, ix]) / d[:, ix], color='r')
    ylim([-0.1, 0.1])
    xlim([0, 1000])
    # xlim([2, 200])


CMBRainbow(r'C:\tmp\Planck\final_Nov14\cl_chains\iso\base_alpha1_plikHM_TT_lowTEB_post', parname='alpha1', cl_ix=3, Lmax=200, ymax=1.3, data=True, tag='EE')
# CMBRainbow(r'C:\tmp\Planck\final_Nov14\cl_chains\iso\base_alpha1_plikHM_TT_lowTEB_post', parname='alpha1', cl_ix=2, Lmax=300, data=True, tag='TE')

CMBRainbow(r'C:\tmp\Planck\final_Nov14\cl_chains\iso\base_alpha1_plikHM_TT_lowTEB_post', parname='alpha1', cl_ix=1, ymax=2000,
            Lmax=60, data=True, tag='TT')

# CMBRainbow(r'C:\tmp\Planck\final_Nov14\cl_chains\base_plikHM_TT_lowTEB_post', parname='ns', cl_ix=3, Lmax=200, data=False)

# CMBRainbow(r'C:\tmp\Planck\final_Nov14\cl_chains\base_nnu_plikHM_TT_lowTEB_post', parname='nnu', cl_ix=2, Lmax=500)

# lensingPaper()

show()
