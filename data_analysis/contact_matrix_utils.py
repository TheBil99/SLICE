import numpy as np
import pickle
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import patches  as patches
from matplotlib import collections  as mc
from scipy.stats import pearsonr, spearmanr
""" Library created by Federico that contains functions useful to analyze GAM matrices"""
""" This library is completely independent from SLICE scripts """

def bp_to_window_conversion(start, stop, resolution):
    inf = round(start/resolution);     
    sup = round(stop/resolution); 
    return inf, sup

def plot_heatmap(mat, vmin = None,  vmax = None, x_figsize = None, y_figsize = None, save_name = None,save = False, dpi = 300, colormap = 'afmhot_r', percent = False, aspect = None, title = None, resolution = None):
    
    if(x_figsize != None and y_figsize != None):
        f, ax = plt.subplots(figsize=(x_figsize, y_figsize))
    else:
        f, ax = plt.subplots(figsize=(7, 7))
    
    if vmin == None or vmax == None:
        #im = plt.imshow(mat, cmap = colormap, interpolation = 'nearest', aspect=aspect)
        im = plt.imshow(mat, cmap = colormap, interpolation = 'nearest', vmin = np.nanpercentile(mat, 3), vmax = np.nanpercentile(mat, 97) , aspect=aspect)
    else:
        im = plt.imshow(mat, cmap = colormap, interpolation = 'nearest', vmin = vmin, vmax = vmax , aspect=aspect )
    c = plt.colorbar(shrink = 0.5)
    if percent == True:
        vals = c.ax.get_yticks();
        c.ax.set_yticklabels(['{:,.2%}'.format(x) for x in vals]);
    if(title != None):
        ax.set_title(title, fontsize = 15)

    if(resolution != None):
        x_ticks = [0, mat.shape[1]]
        x_tick_labels = [0, str(int(np.ceil(mat.shape[1] * resolution / 1E6))) + "Mb"]

        y_ticks = [0, mat.shape[0]]
        y_tick_labels = [0, str(int(np.ceil(mat.shape[0] * resolution / 1E6))) + "Mb"]

        ax.set_xticks(x_ticks);    ax.set_yticks(y_ticks)
        ax.set_xticklabels(x_tick_labels, fontsize = 30);    ax.set_yticklabels(y_tick_labels, fontsize = 30);


    if (save == True and save_name != None):
        plt.savefig(save_name, dpi = dpi, bbox_inches='tight')
    return f, ax, c, im

def plot_genomewide_averages(mat, vmin , vmax , x_figsize= 15, y_figsize= 15, percent=True, colormap='RdYlBu_r', title = None):
    fig, ax, c, _ = plot_heatmap(mat, vmin = vmin, vmax = vmax, x_figsize= x_figsize, y_figsize= y_figsize, percent=percent, colormap = colormap)
    for i in range(19):
        for j in range(19):
            text = ax.text(j, i, str(np.around(mat[i, j]* 100, decimals=1) ) + "%",
                        ha="center", va="center", color="black", fontsize = 13)

    ax.set_title(title, fontsize = 26)
    fig.tight_layout()
    ax.set_xticks(np.arange(0,19, 1))
    ax.set_yticks(np.arange(0,19, 1))
    ax.set_xticklabels(['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19'], fontsize = 14);
    ax.set_yticklabels(['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19'], fontsize = 14);
    c.ax.tick_params(axis='y', which='major', labelsize=20)
    return fig, ax, c

def plot_multiple_heatmaps(n_mats, mats,  vmin = None , vmax = None,x_figsize = None, y_figsize = None,  save_name = None,save = False, dpi = 300, colormap = 'afmhot_r', sharex = False, sharey  = False, put_colorbar = True):


    if(x_figsize != None and y_figsize != None):
        fig, axes = plt.subplots(nrows=1, ncols=n_mats, figsize=(x_figsize, y_figsize), sharex=sharex, sharey = sharey)
    else:
        fig, axes = plt.subplots(nrows=1, ncols=n_mats, figsize=(10, 4), sharex=sharex, sharey = sharey)

    
    heatmap_list = []
    for i in range(n_mats):
        if(vmin != None and vmax !=None):
            heatmap_list.append(axes[i].imshow(mats[i], cmap = colormap, interpolation = 'nearest', vmin = vmin[i], vmax = vmax[i]))
        else:
            heatmap_list.append(axes[i].imshow(mats[i], cmap = colormap, interpolation = 'nearest'))
        if(put_colorbar == True):    
            c = plt.colorbar(heatmap_list[-1], ax = axes[i], shrink = 0.5)
    
    #a = axes[0].imshow(ana_gam, cmap = 'afmhot_r', interpolation = 'nearest')
    #b = axes[1].imshow(in_silico_gam, cmap = 'afmhot_r', interpolation = 'nearest')
    #c = plt.colorbar(a, ax = axes[0], shrink = 0.5)
    #c = plt.colorbar(b, ax = axes[1], shrink = 0.5)

    fig.tight_layout()
    if (save == True and save_name != None):
        plt.savefig(save_name, dpi = dpi, bbox_inches='tight')

    return fig, axes

def plot_heatmap_window_seg_freq(mat, arr_1, arr_2,vmin = None,  vmax = None, x_figsize = 8, y_figsize = 10, save_name = None,save = False, dpi = 300, colormap = 'autumn_r', aspect = 'auto', title = None, resolution = None):
    fig, ax = plt.subplots(2, 2, figsize = (x_figsize,y_figsize), sharex = 'col', sharey = 'row',gridspec_kw={'height_ratios': [ 10,1 ], 'width_ratios' :[10,1]  })

    
    if vmin == None or vmax == None:
        #im = plt.imshow(mat, cmap = colormap, interpolation = 'nearest', aspect=aspect)
        im = ax[0,0].imshow(mat, cmap = colormap, interpolation = 'nearest', vmin = np.nanpercentile(mat, 3), vmax = np.nanpercentile(mat, 97) , aspect=aspect)
    else:
        im = ax[0,0].imshow(mat, cmap = colormap, interpolation = 'nearest', vmin = vmin, vmax = vmax , aspect=aspect )


    ax[1,0].plot(arr_1, color = 'orange')
    ax[1,0].axhline(  np.nanmean(arr_1), color = 'blue', linewidth = 1)

    ax[0,1].plot(arr_2, np.arange(0,arr_2.shape[0], 1), color = 'orange')
    ax[0,1].axvline(  np.nanmean(arr_2), color = 'blue', linewidth = 1)

    ax[0,1].set_xticks([])
    ax[1,1].set_xticks([])
    ax[1,1].set_yticks([])

    ax[1,1].axis('off')

    ax[0,0].set_title(title)

    if(resolution != None):
        x_ticks = [0, mat.shape[1]]
        x_tick_labels = [0, str(int(np.ceil(mat.shape[1] * resolution / 1E6))) + "Mb"]

        y_ticks = [0, mat.shape[0]]
        y_tick_labels = [0, str(int(np.ceil(mat.shape[0] * resolution / 1E6))) + "Mb"]

        ax[1, 0].set_xticks(x_ticks);    ax[0,0].set_yticks(y_ticks)
        ax[1, 0].set_xticklabels(x_tick_labels, fontsize = 30);    ax[0,0].set_yticklabels(y_tick_labels, fontsize = 30);

    plt.subplots_adjust(wspace=0, hspace=0)

    if (save == True and save_name != None):
        plt.savefig(save_name, dpi = dpi, bbox_inches='tight')


    return fig, ax, im

def plot_heatmap_two_chromosomes(mat_chr_A, mat_chr_B, mat_chr_AB, title = None, colormap = "autumn_r", resolution = 1E6):
    matrix_pi_up = np.hstack([mat_chr_A, mat_chr_AB ])
    matrix_pi_down = np.hstack([mat_chr_AB.T,  mat_chr_B])

    matrix_pi = np.vstack([matrix_pi_up, matrix_pi_down])


    #matrix_pi[matrix_pi == 0] = np.nan


    fig, ax, c, _ = plot_heatmap(matrix_pi ,x_figsize=15, y_figsize=15, colormap=colormap);

    chrA_lenght = mat_chr_A.shape[0]
    chrB_lenght = mat_chr_B.shape[0]
    
    ticks = [0, chrA_lenght, chrA_lenght + chrB_lenght]
    ticks_labels = ["0", str(int(np.ceil(chrA_lenght * resolution / 1E6))) + "Mb", str(int(np.ceil(chrB_lenght * resolution / 1E6))) + "Mb" ]





    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks_labels, fontsize = 20)
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticks_labels, fontsize = 20)

    c.ax.tick_params(axis='y', which='major', labelsize=20)

    ax.axvline(chrA_lenght, color = "black");
    ax.axhline(chrA_lenght, color = "black");

    ax.set_title(title, fontsize = 15);


    return fig, ax, c

def draw_square(X, Y, resolution, ax, x_annotation = None,y_annotation = None, linewidth = 1, edgecolor = 'b', linestyle = '-',x_s = 0, y_s = 0):
    """ Function to draw squares on a heatmap, useful to identify selected regions """
    """ The square coordinates have to be in the order x_inf, x_sup, y_inf, y_sup """
    """ x_s and y_s are in case the heatmap does not start from 0 but from another coordinate """

    x_inf = X[0]/resolution - x_s
    x_sup = X[1]/resolution - x_s
    y_inf = Y[0]/resolution - y_s
    y_sup = Y[1]/resolution - y_s

    x_avg = (x_sup + x_inf)/2
    y_avg = (y_sup + y_inf)/2

    rect = patches.Rectangle((x_inf, y_inf), x_sup - x_inf, y_sup - y_inf, linewidth=linewidth, edgecolor=edgecolor, facecolor='none', linestyle = linestyle)
    ax.add_patch(rect)
    # #ax.annotate('', xy=(SE17B[0]/150E3, -1), xytext=(SE17B[1]/150E3, -1), xycoords='data', textcoords='data',arrowprops={'arrowstyle': '-'})
    # ax.annotate(x_annotation, xy=(x_avg, 0), ha='center', va='center', color = 'red')
    # #ax.annotate('', xy=(0,  _14mTAD1[0]/150E3), xytext=(0,  _14mTAD1[1]/150E3), xycoords='data', textcoords='data',  arrowprops={'arrowstyle': '-'})
    # ax.annotate(y_annotation, xy=(0, y_avg), ha='center', va='center', color = 'red')


    # lc = mc.LineCollection([[(x_inf, -1), (x_sup, -1)], [( -1, y_inf), ( -1, y_sup)]  ]  , color='red', linewidths=3)
    # ax.add_collection(lc)
    
    

def add_identity(axes, *line_args, **line_kwargs):
    identity, = axes.plot([], [], *line_args, **line_kwargs)
    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes

def contact_distance_distribution(matrix):
    distribution = np.zeros(matrix.shape[0])
    std          = np.zeros(matrix.shape[0])
    for i in range(matrix.shape[0]):
        temp = np.diag(matrix, i)
        distribution[i] = np.mean(temp)
        std[i]          = np.std(temp)

    return distribution, std

def contact_distance_distribution_nonzero(matrix):
    #tiene conto dei valori a zero e non li inserisce nel calcolo di media e std
    distribution = np.zeros(matrix.shape[0])
    std          = np.zeros(matrix.shape[0])
    for i in range(1, matrix.shape[0]-1):
        temp = np.diag(matrix, i)
        if (temp!=0).sum() != 0:
            distribution[i] = np.true_divide(temp.sum(),(temp!=0).sum())
            std[i]          = (np.true_divide((temp**2).sum(),(temp!=0).sum()) - distribution[i]**2)**(1/2)
        else:
            distribution[i] = 0
            std[i]          = 0
        
        
            
    return distribution, std

def matrix_correlation(mat_1, mat_2, verbose_print = False):
    """ Function that returns pearson and spearman correlation of two matrices containing nan values """


    f_mat   = mat_1.flatten()
    f_mat_1 = mat_2.flatten()
    bad = ~np.logical_or(np.isnan(f_mat), np.isnan(f_mat_1))

    f_mat = np.compress(bad, f_mat)  
    f_mat_1 = np.compress(bad, f_mat_1)  

    p, s = pearsonr(f_mat, f_mat_1)[0], spearmanr(f_mat, f_mat_1)[0]
    
    if verbose_print == True:
        print("pearson\t\t\tspearman\n", p, "\t", s)

    return p,s

def scatter_matrices(mat_1, mat_2, x_label = None, y_label = None, title = None, colors = None):

    fig, ax = plt.subplots(1, figsize = (10,10))
    #if(colors != None):
    ax.scatter(mat_1, mat_2, c = colors);
    # else:
    #     ax.scatter(mat_1, mat_2);
    ax.set_xlabel(x_label, fontsize = 20);
    ax.set_ylabel(y_label, fontsize = 20);
    ax.set_title(title, fontsize = 20);
    ax.tick_params(axis='x', which='major', labelsize=30)
    ax.tick_params(axis='y', which='major', labelsize=30)


    return fig, ax



################################################################
########### deprecated #########################################
################################################################

def select_submatrix(contact_matrix, selection_matrix):
    #scelgo della contact matrix solo le entrate corrispondenti ai valori
    #diversi da zero della selection matrix, per esempio se 
    #voglio selezionare solo i prominent contacts in una 
    #cosegregation matrix
    
    mat = np.copy(selection_matrix)
    mat[np.nonzero(mat)] = mat[np.nonzero(mat)] /mat[np.nonzero(mat)] 
    return contact_matrix*mat

def contact_number_between_selected_windows(contact_matrix, data, dataformat = 'df',  attributes = None):
    #funzione per vedere quanti sono i contatti tra una lista di windows selezionate
    #che possono essere fornite attraverso un dataframe in cui sono individuate le windows di 
    #interesse oppure in una semplice lista
    if(dataformat == 'df'):
        data_temp = data.copy()
        data_temp['selected_windows'] = data_temp.apply(lambda x: x['start_window'] if x['same_window'] == 1 else np.nan, axis = 1)
        selected_windows = data_temp['selected_windows'].values
        selected_windows = selected_windows[~np.isnan(selected_windows)].astype(int)

    if(dataformat == 'list'):
        selected_windows = data

    temp = np.array([len(selected_windows)*[i] for i in selected_windows]).astype(int)
    zero_diagonal_contact_matrix = np.copy(contact_matrix)

    #è importante forzare la diagonale ad essere di zeri se no vengono contati anche i termini in diagonale
    np.fill_diagonal(zero_diagonal_contact_matrix, np.zeros(contact_matrix.shape[0]) )  
    selected_submatrix = zero_diagonal_contact_matrix[temp, selected_windows]

    #diviso due perché i contatti vanno contati solo una volta e la matrice è simmetrica
    
    return np.count_nonzero(selected_submatrix)/2

def contact_number_between_selected_windows_mod(contact_matrix, data, dataformat = 'df',  attributes = None):
    #funzione per vedere quanti sono i contatti tra una lista di windows selezionate
    #che possono essere fornite attraverso un dataframe in cui sono individuate le windows di 
    #interesse oppure in una semplice lista
    if(dataformat == 'df'):
        data_temp = data.copy()
        data_temp['selected_windows'] = data_temp.apply(lambda x: list(range(x['start_window'],  x['end_window']+1)), axis = 1)
        selected_windows = [item for sublist in data_temp['selected_windows'].values for item in sublist]
        

    if(dataformat == 'list'):
        selected_windows = data

    temp = np.array([len(selected_windows)*[i] for i in selected_windows]).astype(int)
    zero_diagonal_contact_matrix = np.copy(contact_matrix)

    #è importante forzare la diagonale ad essere di zeri se no vengono contati anche i termini in diagonale
    np.fill_diagonal(zero_diagonal_contact_matrix, np.zeros(contact_matrix.shape[0]) )  
    selected_submatrix = zero_diagonal_contact_matrix[temp, selected_windows]

    #diviso due perché i contatti vanno contati solo una volta e la matrice è simmetrica
    
    return np.count_nonzero(selected_submatrix)/2

