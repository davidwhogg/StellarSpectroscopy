from scipy import linalg
import numpy as np
import pylab as pl
import matplotlib as mpl
from matplotlib import colors

from sklearn.lda import LDA
from sklearn.qda import QDA

import sqldata_notmain

###############################################################################
# colormap
cmap = colors.LinearSegmentedColormap(
    'red_blue_classes',
    {'red': [(0, 1, 1), (1, 0.7, 0.7)],
     'green': [(0, 0.7, 0.7), (1, 0.7, 0.7)],
     'blue': [(0, 0.7, 0.7), (1, 1, 1)]})
pl.cm.register_cmap(cmap=cmap)


###############################################################################

def dataset_cov():
    array10=sqldata_notmain.sqldata(1)
    magu=array10[0][14]
    magg=array10[0][15]
    magr=array10[0][16]
    magi=array10[0][17]
    maggr=np.array(magg)-np.array(magr)
    magug=np.array(magu)-np.array(magg)
    magri=np.array(magr)-np.array(magi)
    gr=np.array([maggr]).T
    ri=np.array([magri]).T
    ug=np.array([magug]).T
    X=np.hstack((gr,ri)) #pairwise coordinate of array10 values

    m=len(maggr)/2 #SET THESE TWO
    n=len(maggr)-m
    
    y = np.hstack((np.zeros(m), np.ones(n)))
    return X, y


###############################################################################
# plot functions

def plot_data(lda, X, y, y_pred, fig_index):
    splot = pl.subplot(1, 1, fig_index)
    if fig_index == 1:
        pl.title('Linear Discriminant Analysis')
        pl.ylabel('r-i')
        pl.xlabel('g-r')

    tp = (y == y_pred)  # True Positive
    tp0, tp1 = tp[y == 0], tp[y == 1]
    X0, X1 = X[y == 0], X[y == 1]
    X0_tp, X0_fp = X0[tp0], X0[~tp0]
    X1_tp, X1_fp = X1[tp1], X1[~tp1]
    xmin, xmax = X[:, 0].min(), X[:, 0].max()
    ymin, ymax = X[:, 1].min(), X[:, 1].max()


    # class 0: dots
    pl.plot(X0_tp[:, 0], X0_tp[:, 1], 'o', color='red')
    pl.plot(X0_fp[:, 0], X0_fp[:, 1], '.', color='#990000')  # dark red

    # class 1: dots
    pl.plot(X1_tp[:, 0], X1_tp[:, 1], 'o', color='blue')
    pl.plot(X1_fp[:, 0], X1_fp[:, 1], '.', color='#000099')  # dark blue
    
    # class 0 and 1 : areas
    nx, ny = 20, 10
    x_min, x_max = pl.xlim()
    y_min, y_max = pl.ylim()
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, nx),
                         np.linspace(y_min, y_max, ny))
    Z = lda.predict_proba(np.c_[xx.ravel(), yy.ravel()])
    Z = Z[:, 1].reshape(xx.shape)
    pl.pcolormesh(xx, yy, Z, cmap='red_blue_classes', #color gradient background
                  norm=colors.Normalize(0., 1.))
    pl.contour(xx, yy, Z, [0.5], linewidths=2., colors='k') #separation line

    # means ### Extra
    pl.plot(lda.means_[0][0], lda.means_[0][1],
            'o', color='black', markersize=10)
    pl.plot(lda.means_[1][0], lda.means_[1][1],
            'o', color='black', markersize=10)

    return splot
    

## Extra
def plot_ellipse(splot, mean, cov, color):
    v, w = linalg.eigh(cov)
    u = w[0] / linalg.norm(w[0])
    angle = np.arctan(u[1] / u[0])
    angle = 180 * angle / np.pi  # convert to degrees
    # filled gaussian at 2 standard deviation
    ell = mpl.patches.Ellipse(mean, 2 * v[0] ** 0.5, 2 * v[1] ** 0.5,
                              180 + angle, color=color)
    ell.set_clip_box(splot.bbox)
    ell.set_alpha(0.5)
    splot.add_artist(ell)
    splot.set_xticks(())
    splot.set_yticks(())

##Extra
def plot_lda_cov(lda, splot):
    plot_ellipse(splot, lda.means_[0], lda.covariance_, 'red')
    plot_ellipse(splot, lda.means_[1], lda.covariance_, 'blue')



###############################################################################
for i, (X, y) in enumerate([dataset_cov()]):
    # LDA
    lda = LDA()
    y_pred = lda.fit(X, y, store_covariance=True).predict(X)
    splot = plot_data(lda, X, y, y_pred , fig_index=2 * i + 1)
    plot_lda_cov(lda, splot) ##Extra
    #pl.axis('tight')


pl.suptitle('LDA')
pl.show()
