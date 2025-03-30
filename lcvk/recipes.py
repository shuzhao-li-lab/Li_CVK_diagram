#
# Misc. receipes; likely of bad taste
#

def differential_metabolites(list_cpds):
    '''
    Coloring differential metabolites from metabolomics data
    '''
    pass



def plot_global_metabolites(allcpds):
    '''Illustration only'''
    YMAX = 2.5
    theta_limit = (0.1*np.pi, 1.9*np.pi)
    theta_offset = 0.05         # so that nodes do not sit on border lines
    hcList = [max(min(x['ratio_H_C'], 2.5), 0.5) for x in allcpds]
    hcList = standardize_data(hcList)

    # ticks reverse function of standardize_data
    xlabels = [ 0.8, 1.1, 1.4, 1.7, 2, 2.3]
    xticks = standardize_data(xlabels)

    yList = [min(x['ratio_NOPS'] + x['ccs_offset'], YMAX) 
            for x in allcpds]

    fig = plt.figure(figsize=(12, 12))
    ax = fig.subplots(subplot_kw={'projection': 'polar',  } 
    )
    c = ax.scatter(hcList, yList, marker='o', facecolors='w', 
                    linewidths=.3, edgecolors='k', 
                    s=16,
    )
    ax.set_theta_zero_location('S')
    ax.set_theta_direction(-1)
    ax.set_thetalim((theta_limit[0]-theta_offset, theta_limit[1]+theta_offset))
    ax.set_ylim(0, YMAX + 0.1)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xlabels])

    ax.set_rgrids([0.5, 1.5])
    ax.grid(color='r', linewidth=0.1)
    for p in ax.spines.values():
        # inner, polar, start, end
        p.set_color('r')
        p.set_linewidth(0.2)

    ax.set_rorigin(-1)
    ax.text(1.7, YMAX + .45, "H:C", )
    ax.text(1.96*np.pi, 1.6, "NOPS:C", rotation=285)
    ax.set_title("Global metabolites")
    plt.savefig('gem_mar2025_ratio_NOPS.pdf', dpi=300)


def plot_global_metabolites_with_colorExamples(allcpds, 
                        outfile='gem_mar2025_ratio_NOPS_colorExamples.pdf'):
    '''Illustration only'''
    YMAX = 2.5
    theta_limit = (0.1*np.pi, 1.9*np.pi)
    theta_offset = 0.05         # so that nodes do not sit on border lines
    hcList = [max(min(x['ratio_H_C'], 2.5), 0.5) for x in allcpds]
    hcList = standardize_data(hcList)

    # ticks reverse function of standardize_data
    xlabels = [ 0.8, 1.1, 1.4, 1.7, 2, 2.3]
    xticks = standardize_data(xlabels)

    yList = [min(x['ratio_NOPS'] + x['ccs_offset'], YMAX) 
            for x in allcpds]

    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(subplot_kw={'projection': 'polar',  } 
    )

    c = ax.scatter(hcList, yList, marker='o', facecolors='w', 
                    linewidths=.3, edgecolors='k', 
                    s=16,
    )

    collections = [green, magenta, cyan, black]
    colors_collections = ['green', 'magenta', 'cyan', 'blue']
    for ii in range(3):
        XX = [hcList[j] for j in collections[ii]]
        YY = [yList[j] for j in collections[ii]]
        ax.scatter(XX, YY, marker='o', facecolors=colors_collections[ii], s=24, alpha=.5)

    ax.set_theta_zero_location('S')
    ax.set_theta_direction(-1)
    ax.set_thetalim((theta_limit[0]-theta_offset, theta_limit[1]+theta_offset))
    ax.set_ylim(0, YMAX + 0.1)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xlabels])

    ax.set_rgrids([0.5, 1.5])
    ax.grid(color='r', linewidth=0.1)
    for p in ax.spines.values():
        # inner, polar, start, end
        p.set_color('r')
        p.set_linewidth(0.2)

    ax.set_rorigin(-1)
    ax.text(1.7, YMAX + .45, "H:C", )
    ax.text(1.96*np.pi, 1.6, "NOPS:C", rotation=285)
    ax.set_title("Global metabolites")
    plt.savefig(outfile, dpi=300)
    
    
def plot_circular_pathway(pathway_nodes, pathway_edges, title, 
                          ydata='ratio_NOPS', 
                      fontsize=5, rotation=30, show_names=True, show_formula=False, outdir='.'
                      ):
    theta_limit = (0.1*np.pi, 1.9*np.pi)
    theta_offset = 0.05         # so that nodes do not sit on border lines
    hcList = [max(min(x['ratio_H_C'], 2.5), 0.3) for x in pathway_nodes]
    hcList = standardize_data(hcList)
    
    # ticks reverse function of standardize_data
    xlabels = [0.8, 1.1, 1.4, 1.7, 2, 2.3]
    xticks = standardize_data(xlabels)
    
    YMAX = 3
    yList = [min(x[ydata] + x['ccs_offset'], YMAX)  for x in pathway_nodes]
    
    formulas = [x['neutral_formula'] for x in pathway_nodes]
    names = [x['name'] for x in pathway_nodes]
    
    layout = dict((x['id'], (hcList[ii], yList[ii])) for ii,x in enumerate(pathway_nodes))
    
    _nodes = [x['id'] for x in pathway_nodes]
    #sizes = [len(fdict[x[0]]) * 3 for x in kv_list_gem]
    
    text_offset = 0.01
    fig = plt.figure(figsize=(15,15))
    ax = fig.subplots(subplot_kw={'projection': 'polar'})
    ax.set_theta_zero_location('S')
    ax.set_theta_direction(-1)
    ax.set_thetalim((theta_limit[0]-theta_offset, theta_limit[1]+theta_offset))
    ax.set_ylim(0, 3.1)
    
    ax.scatter(hcList, yList, marker='o', facecolors='blue', 
                linewidths=.2, edgecolors='m', 
                s=48, alpha=0.5,
                )
    
    # optionally optimize text labels here, to avoid too much overlap in cpd names
    
    for ii, f in enumerate(formulas):
        L = names[ii]
        if show_names:
            ax.text(hcList[ii], yList[ii] - text_offset, L, color='k',
                 fontsize=fontsize, rotation=rotation, alpha=0.9)
        if show_formula:
            L = f + ', \n' + names[ii]
            ax.text(hcList[ii], yList[ii] - text_offset, L, color='k',
                 fontsize=fontsize, rotation=rotation, alpha=0.9)
        
    for edge in set(pathway_edges):     # set to remove repeated edges; edge has to be tuple
        if edge[0] in _nodes and edge[1] in _nodes:
            ax.annotate("",
                xy=layout[edge[0]], xycoords='data',
                xytext=layout[edge[1]], textcoords='data',
                arrowprops=dict(arrowstyle="->", color="0.5", alpha=0.2,
                                shrinkA=5, shrinkB=5,
                                patchA=None, patchB=None,
                                connectionstyle="arc3,rad=-0.3",
                                ),
                )
    
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xlabels])

    ax.set_rgrids([0.5, 2])
    ax.grid(color='r', linewidth=0.05)
    for p in ax.spines.values():
        # inner, polar, start, end
        p.set_color('r')
        p.set_linewidth(0.1)

    ax.set_rorigin(-1)
    ax.text(1.7, YMAX + .35, "H:C", )
    ax.text(1.945*np.pi, 1.6, ydata, rotation=285)
    
    ax.set_title(title) # "KV-Li plot"
    plt.savefig(
        os.path.join(outdir, 
        title.replace('/', '__') + '.pdf')
        )
    