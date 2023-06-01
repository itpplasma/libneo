#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 12:36:46 2020

@author: ha_pulferer
"""
from numpy import where, linspace, meshgrid, loadtxt
import matplotlib.pyplot as plt
from matplotlib import cm


def interpol(solution):
    """
    Interpolate bounce or drift frequency over constants of motion (pph, H)
    for nr random initial values
    """
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.interpolate import Rbf
    
    x = solution[:,2]       # pphi
    y = solution[:,3]       # H
    z = solution[:,4]       # mu
    w = solution[:,0]       # w_bounce
    type_ = solution[:,5]   # orbit type (passing=0/trapped=1)
    idx_pas = where(type_ == 0)[0]
    idx_ban = where(type_ == 1)[0]
    idx_kid = where(type_ == 2)[0]
    idx_ck = where(type_ == 3)[0]
    idx_oc = where(type_ == 4)[0]
    idx_ic = where(type_ == 5)[0]
    
    # New values for interpolation
    yi = linspace(min(y),max(y),200)
    zi = linspace(min(z),max(z),200)
    Yi,Zi = meshgrid(yi,zi)
    
    # 2D Rbf interpolation in pphi and mu
    rbf_lin = Rbf(y,z,w,function='linear')
    W_rbf_lin = rbf_lin(Yi,Zi)
    
    
    ### Plotting ###
    
    # 2D plot (pphi, mu, w_bounce)
    plt.ticklabel_format(style='sci', scilimits=(0,0))
    plt.title('Radial Basis Functions lin for $\omega_b$')
    plt.pcolor(Yi,Zi,W_rbf_lin,cmap=cm.jet)
    plt.xlabel('H')
    plt.ylabel('$\mu$')
    plt.colorbar()
    plt.show()
    
    # 3D plot (pphi,mu,w_bounce)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('$p_\phi$')
    ax.set_ylabel('$\mu$')
    ax.set_zlabel('$\omega_b$')
    ax.ticklabel_format(style='sci', scilimits=(0,0))
    if len(idx_pas) != 0:
        ax.scatter(x[idx_pas], z[idx_pas], w[idx_pas], marker='.', label='Passing Orbits')
    elif len(idx_ban) != 0:
        ax.scatter(x[idx_ban], z[idx_ban], w[idx_ban], marker='v', label='Banana Orbits')
    elif len(idx_ban) != 0:
        ax.scatter(x[idx_ban], z[idx_ban], w[idx_ban], marker='*', label='Potato Orbits')
    elif len(idx_ban) != 0:
        ax.scatter(x[idx_ban], z[idx_ban], w[idx_ban], marker='+', label='Concave-Kidney Orbits')
    elif len(idx_ban) != 0:
        ax.scatter(x[idx_ban], z[idx_ban], w[idx_ban], marker='x', label='Outer-Circulating Orbits')
    elif len(idx_ban) != 0:
        ax.scatter(x[idx_ban], z[idx_ban], w[idx_ban], marker='D', label='Inner-Circulating Orbits')
    ax.legend()
    plt.show()

    
    # 2D plot - heatmap (pphi, mu, w_bounce (trapped/passing))
    for axs in enumerate([x,y]):
        fig, ax = plt.subplots()
        if len(idx_pas) != 0:
            img21 = ax.scatter(axs[1][idx_pas], z[idx_pas], c=w[idx_pas], marker='.', cmap='Reds', label='Passing Orbits')
            cbar21 = fig.colorbar(img21)
            cbar21.set_label('Bounce Frequency $\omega_b$ / Hz')
        if len(idx_ban) != 0:
            ax.scatter(axs[1][idx_ban], z[idx_ban], c=w[idx_ban], marker='v', cmap='Reds', label='Banana Orbits')
        if len(idx_kid) != 0:
            ax.scatter(axs[1][idx_kid], z[idx_kid], c=w[idx_kid], marker='*', cmap='Reds', label='Potato Orbits')
        if len(idx_ck) != 0:
            ax.scatter(axs[1][idx_ck], z[idx_ck], c=w[idx_ck], marker='+', cmap='Reds', label='Concave-Kidney Orbits')
        if len(idx_oc) != 0:
            ax.scatter(axs[1][idx_oc], z[idx_oc], c=w[idx_oc], marker='x', cmap='Reds', label='Outer-Circulating Orbits')
        if len(idx_ic) != 0:
            ax.scatter(axs[1][idx_ic], z[idx_ic], c=w[idx_ic], marker='D', cmap='Reds', label='Inner-Circulating Orbits')
        plt.xlim([min(axs[1]) - (max(axs[1])-min(axs[1]))*0.1, max(axs[1]) + (max(axs[1])-min(axs[1]))*0.1])
        plt.ylim([min(z) - (max(z)-min(z))*0.1, max(z) + (max(z)-min(z))*0.1])
        if axs[0] == 0: plt.xlabel('$p_\phi$')
        else: plt.xlabel('H')
        plt.ylabel('$\mu$')
        ax.legend()
        plt.show()
     
    for axs in enumerate([x,y,z]):
        fig, ax = plt.subplots()
        if len(idx_pas) != 0:
            plt.plot(axs[1][idx_pas], w[idx_pas], 'r.', label='Passing Orbits')
        if len(idx_ban) != 0:
            plt.plot(axs[1][idx_ban], w[idx_ban], 'rv', label='Banana Orbits')
        if len(idx_kid) != 0:
            plt.plot(axs[1][idx_kid], w[idx_kid], 'r*', label='Potato Orbits')
        if len(idx_ck) != 0:
            plt.plot(axs[1][idx_ck], w[idx_ck], 'r+', label='Concave-Kidney Orbits')
        if len(idx_oc) != 0:
            plt.plot(axs[1][idx_oc], w[idx_oc], 'rx', label='Outer-Circulating Orbits')
        if len(idx_ic) != 0:
            plt.plot(axs[1][idx_ic], w[idx_ic], 'rD', label='Inner-Circulating Orbits')
        plt.ylabel('$\omega_b$ / Hz')
        if axs[0] == 0: plt.xlabel('$p_{\phi}$')
        elif axs[0] == 1: plt.xlabel('H')
        elif axs[0] == 2: plt.xlabel('$\mu$')
        plt.legend()
        fig.show()
        
    

def suppvec(solution):
    
    import numpy as np
    from sklearn.model_selection import train_test_split    
    from sklearn import svm
    
    X = solution[:,[2,4]]*1e35   # wb, pphi, H, mu
    y = solution[:,5]           # orbit type
    y = y.astype(int)
    
    def make_meshgrid(x,y,nr):
        x_min, x_max = x.min() - 0.1*(x.max()-x.min()), x.max() + 0.1*(x.max()-x.min())
        y_min, y_max = y.min() - 0.1*(y.max()-y.min()), y.max() + 0.1*(y.max()-y.min())
        xx, yy = np.meshgrid(np.linspace(x_min,x_max,nr), np.linspace(y_min,y_max,nr))
        return xx, yy
    
    def plot_contours(ax, clf, xx, yy, **params):
        Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)
        out = ax.contourf(xx,yy,Z,**params)
        return out
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.4)

    clf = svm.SVC(kernel='linear',decision_function_shape='ovo', C=100)
    clf.fit(X_train, y_train)
    

    fig, ax = plt.subplots()
    
    X0, X1 = X[:,0], X[:,1]
    xx, yy = make_meshgrid(X0,X1,100)
    
    plot_contours(ax, clf, xx, yy, cmap=plt.cm.coolwarm, alpha=0.8)
    ax.scatter(X0,X1,c=y,cmap=plt.cm.coolwarm,s=20,edgecolors='k')
    ax.set_xlabel('$\omega_b$')
    ax.set_ylabel('$p_\phi$')
    ax.legend()
    plt.show()
    
def suppvec2(solution):
    
    from matplotlib.colors import ListedColormap
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn import svm
    from sklearn.model_selection import train_test_split
    
    def plot_decision_regions(X, y, classifier, nr = 200):
    
        # setup marker generator and color map
        markers = ('s', '.')#, '+', '^','v')#'o', '^', 'x', '+', 'v')
        colors = ('darkorange', 'orangered')#, 'red', 'teal','aqua') #'orange', 'purple', 'red', 'black')
        labels = ('passing', 'banana')#, 'kidney', 'concave-kidney', 'outer-circulating') #, 'inner-circulating')
        cmap = ListedColormap(colors[:len(np.unique(y))+1])
    
        # plot the decision surface
        x1_min, x1_max = X[:, 0].min() - 0.1*( X[:, 0].max() - X[:, 0].min() ), X[:, 0].max() + 0.1*( X[:, 0].max() - X[:, 0].min() )
        x2_min, x2_max = X[:, 1].min() - 0.1*( X[:, 1].max() - X[:, 1].min() ), X[:, 1].max() + 0.1*( X[:, 1].max() - X[:, 1].min() )
        xx1, xx2 = np.meshgrid(np.linspace(x1_min, x1_max, nr),
                               np.linspace(x2_min, x2_max, nr))
        Z = classifier.predict(np.array([xx1.ravel(), xx2.ravel()]).T)
        Z = Z.reshape(xx1.shape)
        plt.contourf(xx1, xx2, Z, alpha=0.4, cmap=cmap)
        plt.xlim(xx1.min(), xx1.max())
        plt.ylim(xx2.min(), xx2.max())
    
        # plot class samples
        for idx, cl in enumerate(np.unique(y)):
            plt.scatter(x=X[y == cl, 0], y=X[y == cl, 1],
                        alpha=0.8, c=colors[idx],
                        marker=markers[idx], label=labels[idx])
            plt.xlabel('$\\bar{r}$ / cm')
            plt.ylabel('$\lambda_0 = \mu B_0 / H$')
            plt.legend()
            plt.title('Regions for different orbit types (v = $9.8e5$ m/s)')
    
    
    # Loading some example data
#    iris = datasets.load_iris()
#    X = iris.data[:, [0,2]]
    X = solution[:,[7,6]]
#    y = iris.target
    y = solution[:,5].astype(int)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)
#    y = np.concatenate((y, np.ones(50)+2))
#    X = np.concatenate((X, X[:50]*2))
    
#    lr = LogisticRegression(solver='newton-cg', multi_class='multinomial')
    lr =  svm.SVC(C=100, decision_function_shape='ovo')
    lr.fit(X_train, y_train)
    plot_decision_regions(X, y, classifier=lr)
    
    
     
    
    
    
""" ################################## Try ################################ """

solution = loadtxt('solution_2306.csv', delimiter=',')

#interpol(solution)
#plt.show()
suppvec2(solution)

















