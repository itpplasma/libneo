from pylab import arange,pi,sin,cos,sqrt
import matplotlib.pylab as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from matplotlib.ticker import ScalarFormatter,MaxNLocator
from math import log10, floor
from pylab import get_cmap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import TransformedBbox, Affine2D
from IPython import embed


PAD_INCHES = 0.1


def tight_layout(pad_inches=PAD_INCHES, h_pad_inches=None, w_pad_inches=None):
    """Adjust subplot parameters to give specified padding.

    Parameters
    ----------
    pad_inches : float
        minimum padding between the figure edge and the edges of subplots.
    h_pad_inches, w_pad_inches : float
        minimum padding (height/width) between edges of adjacent subplots.
        Defaults to `pad_inches`.
    """
    if h_pad_inches is None:
        h_pad_inches = pad_inches
    if w_pad_inches is None:
        w_pad_inches = pad_inches
    fig = plt.gcf()
    tight_borders(fig, pad_inches=pad_inches)
    # NOTE: border padding affects subplot spacing; tighten border first
    tight_subplot_spacing(fig, h_pad_inches, w_pad_inches)


def tight_borders(fig, pad_inches=PAD_INCHES):
    """Stretch subplot boundaries to figure edges plus padding."""
    # call draw to update the renderer and get accurate bboxes.
    fig.canvas.draw()
    bbox_original = fig.bbox_inches
    bbox_tight = _get_tightbbox(fig, pad_inches)

    # figure dimensions ordered like bbox.extents: x0, y0, x1, y1
    lengths = np.array([bbox_original.width, bbox_original.height,
                        bbox_original.width, bbox_original.height])
    whitespace = (bbox_tight.extents - bbox_original.extents) / lengths

    # border padding ordered like bbox.extents: x0, y0, x1, y1
    current_borders = np.array([fig.subplotpars.left, fig.subplotpars.bottom,
                                fig.subplotpars.right, fig.subplotpars.top])

    left, bottom, right, top = current_borders - whitespace
    fig.subplots_adjust(bottom=bottom, top=top, left=left, right=right)


def _get_tightbbox(fig, pad_inches):
    renderer = fig.canvas.get_renderer()
    bbox_inches = fig.get_tightbbox(renderer)
    return bbox_inches.padded(pad_inches)


def tight_subplot_spacing(fig, h_pad_inches, w_pad_inches):
    """Stretch subplots so adjacent subplots are separated by given padding."""
    # Zero hspace and wspace to make it easier to calculate the spacing.
    fig.subplots_adjust(hspace=0, wspace=0)
    fig.canvas.draw()

    figbox = fig.bbox_inches
    ax_bottom, ax_top, ax_left, ax_right = _get_grid_boundaries(fig)
    nrows, ncols = ax_bottom.shape

    subplots_height = fig.subplotpars.top - fig.subplotpars.bottom
    if nrows > 1:
        h_overlap_inches = ax_top[1:] - ax_bottom[:-1]
        hspace_inches = h_overlap_inches.max() + h_pad_inches
        hspace_fig_frac = hspace_inches / figbox.height
        hspace = _fig_frac_to_cell_frac(hspace_fig_frac, subplots_height, nrows)
        fig.subplots_adjust(hspace=hspace)

    subplots_width = fig.subplotpars.right - fig.subplotpars.left
    if ncols > 1:
        w_overlap_inches = ax_right[:,:-1] - ax_left[:,1:]
        wspace_inches = w_overlap_inches.max() + w_pad_inches
        wspace_fig_frac = wspace_inches / figbox.width
        wspace = _fig_frac_to_cell_frac(wspace_fig_frac, subplots_width, ncols)
        fig.subplots_adjust(wspace=wspace)


def _get_grid_boundaries(fig):
    """Return grid boundaries for bboxes of subplots

    Returns
    -------
    ax_bottom, ax_top, ax_left, ax_right : array
        bbox cell-boundaries of subplot grid. If a subplot spans cells, the grid
        boundaries cutting through that subplot will be masked.
    """
    nrows, ncols, n = fig.axes[0].get_geometry()
    # Initialize boundaries as masked arrays; in the future, support subplots
    # that span multiple rows/columns, which would have masked values for grid
    # boundaries that cut through the subplot.
    ax_bottom, ax_top, ax_left, ax_right = [np.ma.masked_all((nrows, ncols))
                                            for n in range(4)]
    renderer = fig.canvas.get_renderer()
    px2inches_trans = Affine2D().scale(1./fig.dpi)
    for ax in fig.axes:
        ax_bbox = ax.get_tightbbox(renderer)
        x0, y0, x1, y1 = TransformedBbox(ax_bbox, px2inches_trans).extents
        nrows, ncols, n = ax.get_geometry()
        # subplot number starts at 1, matrix index starts at 0
        i = n - 1
        ax_bottom.flat[i] = y0
        ax_top.flat[i] = y1
        ax_left.flat[i] = x0
        ax_right.flat[i] = x1
    return ax_bottom, ax_top, ax_left, ax_right


def _fig_frac_to_cell_frac(fig_frac, subplots_frac, num_cells):
    """Return fraction of cell (row/column) from a given fraction of the figure

    Parameters
    ----------
    fig_frac : float
        length given as a fraction of figure height or width
    subplots_frac : float
        fraction of figure (height or width) occupied by subplots
    num_cells : int
        number of rows or columns.
    """
    # This function is reverse engineered from the calculation of `sepH` and
    # `sepW` in  `GridSpecBase.get_grid_positions`.
    return (fig_frac * num_cells) / (subplots_frac - fig_frac*(num_cells-1))

def multipleSubplot():
    params = {\
                'axes.labelsize': 12,\
                'text.fontsize': 12,\
                'axes.titlesize': 12,\
                'legend.fontsize': 12,\
                'xtick.labelsize': 10,\
                'ytick.labelsize': 10,\
                'pdf.fonttype' : 42,
                 }
    plt.rcParams.update(params)

def fullScreen():
    mng = plt.get_current_fig_manager()
    new_size = (mng.window.maxsize()[0],mng.window.maxsize()[1]-32)
    mng.resize(*new_size)
    #mng.frame.Maximize(True)

def usetex(height=0,width=372.):
  #width=600 for one page landscape figures
  fig_width_pt=width# Get this from LaTeX using \showthe\columnwidth
  inches_per_pt = 1.0/72.27               # Convert pt to inch
  golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
  fig_width = fig_width_pt*inches_per_pt  # width in inches
  if height==0:
    fig_height=fig_width_pt*inches_per_pt*golden_mean    # height in inches
  else:
    fig_height=height*inches_per_pt
  fig_size =  [fig_width,fig_height]
  #fig_size=[2.5736820257368205,1.5906229681400368]#half size
  #fig_size=[3.4315760343157606,2.120830624186716] #2/3

  print 'fig size:', fig_size
  params = {#'backend': 'pdf',
            'axes.labelsize': 14,
            'text.fontsize': 14,
            'axes.titlesize': 14,
            'legend.fontsize': 8,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'figure.figsize': fig_size,
            'axes.grid': False,
            #'text.usetex': True,
            'legend.frameon':False,
            'legend.borderpad': 0.2,
            'legend.labelspacing': 0.2,
            'legend.handletextpad': 0.2,
            'legend.fancybox':False,
            'mathtext.default':'sf',
            'tick.direction':'in',
            'font.family':'sans',
            'font.sans':['Helvetica'],
            #'font.monospace':['Andale Mono'],
            #'font.variant':'it',
            #'lines.antialiased':'True',
            'lines.linewidth':1.0,
            'patch.linewidth':1.0,
            'axes.linewidth':1,
            'xtick.major.width':1,
            'ytick.major.width':1,
            'tick.linewidth':1,
            'tick.length':1,
            'tick.scilimits':(1,-1),
            'xtick.major.pad':5,
            'xtick.minor.pad':5,
            'ytick.major.pad':5,
            'ytick.minor.pad':5,
            'mathtext.fontset':'cm',
            'pdf.fonttype' : 42,
            }
  plt.rcParams.update(params)

def useslide(height=0,width=372.):
  #latex scale min=0.5
  #width=307.28987 beamer
  #width=600 for one page landscape figures
  fig_width_pt=width# Get this from LaTeX using \showthe\columnwidth
  inches_per_pt = 1.0/72.27               # Convert pt to inch
  golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
  fig_width = fig_width_pt*inches_per_pt  # width in inches
  if height==0:
    fig_height=fig_width_pt*inches_per_pt*golden_mean    # height in inches
  else:
    fig_height=height*inches_per_pt
  fig_size =  [fig_width,fig_height]
  #fig_size=[2.5736820257368205,1.5906229681400368]#half size
  #fig_size=[3.4315760343157606,2.120830624186716] #2/3

  params = {#'backend': 'pdf',
            'axes.labelsize': 26,
            'text.fontsize': 26,
            'axes.titlesize': 26,
            'legend.fontsize': 18,
            'xtick.labelsize': 18,
            'ytick.labelsize': 18,
            'figure.figsize': fig_size,
            'axes.grid': False,
            #'text.usetex': True,
            'legend.frameon':False,
            'legend.borderpad': 0.5,
            'legend.labelspacing': 0.1,
            'legend.handletextpad': 0.1,
            'legend.fancybox':False,
            'mathtext.default':'sf',
            'tick.direction':'in',
            'font.family':'sans',
            'font.sans':['Helvetica'],
            #'font.monospace':['Andale Mono'],
            #'font.variant':'it',
            #'lines.antialiased':'True',
            'lines.linewidth':1.5,
            'patch.linewidth':1.5,
            'axes.linewidth':1.5,
            'xtick.markersize':1.5,
            'xtick.major.width':1.5,
            'ytick.major.width':1.5,
            'tick.width':1.0,
            'xtick.length':2.0,
            'tick.scilimits':(3,-2),
            'xtick.major.pad':12,
            'xtick.minor.pad':12,
            'ytick.major.pad':12,
            'ytick.minor.pad':12,
            'mathtext.fontset':'cm',
            'pdf.fonttype' : 42,
            }
  plt.rcParams.update(params)

def on_click(event):
    """
    Enlarge or restore the selected axis.
    Usage: fig.canvas.mpl_connect('button_press_event', fconf.on_click)
    """
    ax = event.inaxes
    full_screen_pos = (0.1, 0.1, 0.85, 0.85)
    #print ax.get_position().bounds
    #if (np.array(full_screen_pos) - np.array(ax.get_position().bounds)) <= 0.01:
        #print 'Works'
    if ax is None:
        # Occurs when a region not in an axis is clicked...
        return
    if event.button is 3:
        # On left click, zoom the selected axes
        if not (np.abs(np.array(full_screen_pos) - np.array(ax.get_position().bounds))).sum() <= 0.01:
            ax._orig_position = ax.get_position()
            ax.set_position([0.1, 0.1, 0.85, 0.85])
            for axis in event.canvas.figure.axes:
            # Hide all the other axes...
                if axis is not ax:
                    axis.set_visible(False)
    #elif event.button is 3:
        # On right click, restore the axes
        else:
            try:
                ax.set_position(ax._orig_position)
                for axis in event.canvas.figure.axes:
                    axis.set_visible(True)
            except AttributeError:
            ## If we haven't zoomed, ignore...
                pass
    else:
        # No need to re-draw the canvas if it's not a left or right click
        return
    event.canvas.draw()

def colors(numcolors,map='spectral'):
    std_col = ['r','b','g','m']
    if numcolors <= 4:
        return std_col[:numcolors]
    cm=plt.get_cmap(map)
    col = []
    for i in range(numcolors):
        col.append(cm(1.*i/numcolors))
    return col

def shapes(nshapes):
    std_shapes = ['o','v','s','*','>','x','d',\
            'p','h','H','D','+','|','_']
    return std_shapes[:nshapes]

def plt_vessel(ax, pol=True, polygon=False):
#    
    if polygon:
        from matplotlib.patches import Polygon

    #embed()
    if pol:
        f = open("/afs/ipp/u/mwillens/repository/dat/pol_vessel.data",'r')
        print 'pol_vessel.data'
    else:
        f = open("/afs/ipp/u/mwillens/repository/dat/tor_vessel.data",'r')

    lines = f.readlines()
    f.close()
    
    vessel = []
    r=[]
    z=[]
    for line in lines:
        if '#' in line:
            continue
        if line == '\n':
            if len(r) > 1:
                if pol:
                    vessel.append({'r':r,'z':z})
                else:
                    vessel.append({'r':np.array(r)*np.cos(np.deg2rad(-22.5*3.))-np.array(z)*np.sin(np.deg2rad(-22.5*3.)),'z':np.array(r)*np.sin(np.deg2rad(-22.5*3.))+np.array(z)*np.cos(np.deg2rad(-22.5*3.)) })
            r=[]
            z=[]
            continue
        #embed()
#try:
        val = line.split()
        if len(val)>1:
            r.append(float(val[0]))
            z.append(float(val[1]))
      
            #except:
            #print 'uppsi'
            #embed()
    #embed()
    

    #if pol==False:
        
        #rotate to new coordinate system
    #    xnew=(np.array(r)*np.cos(np.deg2rad(-22.5*3.))-np.array(z)*np.sin(np.deg2rad(-22.5*3.))).tolist()
    #    ynew = (np.array(r)*np.sin(np.deg2rad(-22.5*3.))+np.array(z)*np.cos(np.deg2rad(-22.5*3.))).tolist()
        
    #    embed()
    #    r=xnew
    #    z=ynew

    for key in range(len(vessel)):
        if polygon:
            ax.add_patch(Polygon(zip(vessel[key]['r'],vessel[key]['z']),facecolor='grey',edgecolor='none'))
        else:        
            ax.plot(vessel[key]['r'],vessel[key]['z'],'k-')
        


def plt_eq_pol(ax,nshot,time,rhos=np.linspace(0.,1.,num=10),linewidth=1,linecolor='k',linestyleSep='solid',linestyleFlux='dashed',diag='EQI',exp='AUGD',ed=0):
    import kk
    kk = kk.KK()
    out1 = kk.kkrhopto(nshot,time,rhos,diag=diag,exp=exp,ed=ed)
    
    jrho = 0
    if np.size(np.squeeze(out1.pf))>1:
        for pf in np.squeeze(out1.pf):
            out2 = kk.kkeqpsp(nshot,time,pf,diag=diag,exp=exp,ed=ed)
            if rhos[jrho] == 1.:
                ax.plot(out2.r_surf,out2.z_surf,linestyle=linestyleSep,color=linecolor,linewidth=linewidth)
            else:
                ax.plot(out2.r_surf,out2.z_surf,linestyle=linestyleFlux,color=linecolor,linewidth=linewidth)
            jrho += 1
    else:
        out2 = kk.kkeqpsp(nshot,time,out1.pf,diag=diag,exp=exp,ed=ed)
       
        if out2.err == 0:
            ax.plot(out2.r_surf,out2.z_surf,linestyle=linestyleFlux,color=linecolor,linewidth=linewidth)
        else:
            r=[]
            z=[]
            for a in np.linspace(0,360,360):
                out3 = kk.kkrhorz(nshot,time,rhos,angle=float(a),diag=diag,exp=exp,ed=ed)
                r.append(out3.r)
                z.append(out3.z)

            ax.plot(r,z,linestyle=linestyleFlux,color=linecolor,linewidth=linewidth)

def plt_eq_tor(ax,nshot,time,rhos=np.linspace(0.,1.,num=10),linewidth=1,linecolor='k',linestyleSep='solid',linestyleFlux='dashed',diag='EQI',exp='AUGD'):
    from matplotlib.patches import Circle
    import kk
    kk = kk.KK()
    out1 = kk.kkrhopto(nshot,time,rhos,diag=diag,exp=exp)
    jrho = 0
    if np.size(out1.pf):
        out2 = kk.kkeqpsp(nshot,time,out1.pf,diag=diag,exp=exp)
        ax.add_patch(Circle((0.,0.),np.max(out2.r_surf),fill=False,edgecolor='k'))
    else:
        for pf in out1.pf:
            out2 = kk.kkeqpsp(nshot,time,pf,diag=diag,exp=exp)
            if rhos[jrho] == 1.:
                ax.add_patch(Circle((0.,0.),np.max(out2.r_surf),fill=False,edgecolor='k'))
            else:
                ax.add_patch(Circle((0.,0.),np.max(out2.r_surf),linestyle='dashed',fill=False,edgecolor='k'))
            jrho += 1




#rm -r ~/.matplotlib/*cache
