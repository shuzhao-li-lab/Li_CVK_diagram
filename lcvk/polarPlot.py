# Functions for Li-CVK pathway diagrams
import matplotlib.pyplot as plt
import numpy as np
from mass2chem.formula import parse_chemformula_dict

# used for optimization
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from scipy.spatial.distance import euclidean

# import networkx as nx

def get_chnops_ratios(f):
    '''
    f: formula
    Returns H/C, O/C, N/C and weighted NOPS/C ratios
    '''
    _d = parse_chemformula_dict(f)
    if 'C' in _d:
        hc = _d['H']/_d['C']
        nc, oc, pc, sc  = _d.get('N', 0), _d.get('O', 0), _d.get('P', 0), _d.get('S', 0)
        nops_c = (nc*14 + oc*16 + pc*31 + sc*32) /(_d['C'] * 12)
        oc, nc = oc/_d['C'], nc/_d['C']
        return hc, oc, nc, nops_c
    else:
        return [0] * 4


def get_isotopic_chnops_ratios(f):
    '''To implement'''
    pass



def get_kvli_ratios(f, ccs_offset):
    '''
    Returns H/C, NOPS coefficient, (weighted NOPS/C), CCS offset (CCS deviation from all isomers in pathway).
    ccs_deviation to be calculated elsewhere, as delta to average CCS of all isomers in model.
    '''
    hc, oc, nc, nops_c = get_chnops_ratios(f)
    return hc, oc, nc, nops_c, ccs_offset

def project_hc_data_radial(hcList, 
                           min_theta=0.1*np.pi, max_theta=1.9*np.pi,
                           ):
    '''
    Convert H:C ratios to radial theta angles.
    min_value=0.5, max_value=2.5, 
    '''
    _low, _high = min(hcList), max(hcList)
    range = _high - _low
    if range==0:
        return hcList
    else:
        _slope = (max_theta - min_theta)/range
        return [_slope*(x-_low) + min_theta for x in hcList]

def standardize_data(list_hc_data, range_degrees=1.8):
    '''
    Use fixed scale. Use at risk. 
    Theta ranges 0.1*pi to 1.9*pi. H:C btw 0.5 and 2.5.
    Therefore HC = 2.5703939893 * x + 0.3
    '''
    return [2.827433388230814 * (x-0.5) + 0.314159265359 for x in list_hc_data]

def clean_pathway_name(name):
    return name.replace('/', '_')

def calculate_clean_network(list_cpds, list_edges, formula='neutral_formula', max_H_C_allowed=5):
    pathway_nodes, pathway_edges = [], []
    for cpd in list_cpds:
        if formula in cpd and 'C' in cpd[formula] and 'H' in cpd[formula]:
            hc, oc, nc, nops_c, ccs_offset = get_kvli_ratios(cpd[formula], 0)
            if hc <= max_H_C_allowed:
                cpd.update(
                            {
                                'ratio_H_C': hc,
                                'ratio_O_C': oc,
                                'ratio_N_C': nc,
                                'ratio_NOPS': nops_c,
                                'ccs_offset': ccs_offset,
                                'centrality': None, # future use
                            }
                        )
                pathway_nodes.append(cpd)
    valid_cpd_ids = [x['id'] for x in pathway_nodes]
    for edge in list_edges:
        if edge[0] in valid_cpd_ids and edge[1] in valid_cpd_ids:
            pathway_edges.append(tuple(edge))
            
    return pathway_nodes, set(pathway_edges)


def get_pathway_coordinates(pathway, cpdDict, formula='charged_formula'):
    '''
    Input a pathway in mummichog style.
    formula : neutral_formula, charged_formula
    Return coordinates {cpd: (hc, nops_c, ccs_offset), ...}
    
    Maybe add ccs_offset later.
    '''
    cpds = []
    for C in pathway['cpds']:
        if C in cpdDict and cpdDict[C][formula]:
            try:
                f = cpdDict[C][formula]
                if 'C' in f and 'H' in f:
                    hc, oc, nc, nops_c, ccs_offset = get_kvli_ratios(f, 0) # not using ccs_offset for now
                    cpds.append(
                        {
                            'id': cpdDict[C]['id'],
                            'name': cpdDict[C]['name'],
                            formula: f,
                            'ratio_H_C': hc,
                            'ratio_O_C': oc,
                            'ratio_N_C': nc,
                            'ratio_NOPS': nops_c,
                            'ccs_offset': ccs_offset,
                            'centrality': None, # future use
                        }
                    )
                else:
                    print(C, f)
            except KeyError:
                print(C)
    return cpds

def generate_landmarks(hcList, yList):
    '''returns the geo center and the max of input coordinates'''
    # Not used now, since there's usually conversion btw Cartesian and Polar
    return (np.mean(hcList), np.mean(yList)), (max(hcList), max(yList))

#
# Main function for Li-cvK diagram
#
def cplot_LCVK_pathway(list_cpds, list_edges, 
                    formula='charged_formula', cpd_name='name',
                      ydata='ratio_NOPS', 
                      show_names=True, show_formula=False, 
                      flexible_theta=True,
                      show_landmarks=False, 
                      max_H_C_allowed=5, 
                      min_theta=0.1*np.pi, max_theta=1.9*np.pi,
                      min_radius=0, max_radius=2.5,
                      padding=0.1, optimization=False,
                      optimization_spacer=0.3, optimization_iter=3,
                      fontsize=5, rotation=30, 
                      width=15, height=15,   # inch
                      yLabel_off=False,
                      # marker parameters 
                      marker='o', facecolors='blue', linewidths=.2, edgecolors='m', s=48, alpha=0.5,
                      title='', outfile='lcvk_plot.pdf', dpi=300
                      ):
    '''
    formula='charged_formula', but recommend neutral
    cpd_name can use 'id'
    '''
    pathway_nodes, pathway_edges = calculate_clean_network(list_cpds, list_edges, formula, 
                                                           max_H_C_allowed)
    hcList = [x['ratio_H_C'] for x in pathway_nodes]
    cartesian_mean_x, cartesian_min_x, cartesian_max_x = np.mean(hcList), min(hcList), max(1, max(hcList))
    hcList = project_hc_data_radial(hcList)
    if flexible_theta:
        min_theta, max_theta = min(hcList), max(hcList)
    # tick labels; showing original data before polar projection 
    xlabels = [ 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4]
    xlabels = [x for x in xlabels if cartesian_min_x < x < cartesian_max_x]
    if not xlabels:  # data falling out of range, forcing 
        xlabels = [1.0, 1.1]
    xticks = project_hc_data_radial(xlabels)
    
    ylabels = {'ratio_NOPS': 'NOPS:C', 'ratio_O_C': 'O:C', 'ratio_H_C': 'H:C'} # proper text for axis labels
    yList = [min(max(x[ydata], min_radius), max_radius) 
             for x in pathway_nodes]
    formulas = [x[formula] for x in pathway_nodes]
    names = [x[cpd_name] for x in pathway_nodes]
    
    # optionally optimize text labels here, to avoid too much overlap in cpd names
    if optimization:
        hcList, yList = optimize_list_coordinates(hcList, yList, 
                                              distance_cut=0.1, xspacer=optimization_spacer, 
                                              rspacer=optimization_spacer, 
                                              optimization_iter=optimization_iter,
                                              )
    
    layout = dict((x['id'], (hcList[ii], yList[ii])) for ii,x in enumerate(pathway_nodes))
    
    fig = plt.figure(figsize=(width, height))
    ax = fig.subplots(subplot_kw={'projection': 'polar'})
    ax.set_theta_zero_location('S')
    ax.set_theta_direction(-1)
    ax.set_thetalim((min_theta-padding, max_theta+padding))
    
    ax.scatter(hcList, yList, marker=marker, facecolors=facecolors, 
                linewidths=linewidths, edgecolors=edgecolors, s=s, alpha=alpha,
                )
    
    if show_landmarks: # show landmarks if desired
        center_x, center_r, min_x, max_r = np.mean(hcList), np.mean(yList), min(hcList), max(yList)
        ax.scatter([center_x, min_x], [center_r, max_r], marker='X', facecolors='w', 
                    linewidths=.5, edgecolors='r', s=128)
        # labels need to be original H:C values
        ax.text(5.7, max_r+.5, f"H:C marks\n{cartesian_min_x:.2f}, {cartesian_mean_x:.2f}", color='r', alpha=.5)
        # ax.text(center_x, center_r-padding, f"({cartesian_mean_x:.2f}, {center_r:.2f})", color='r', alpha=0.5)
        # ax.text(min_x, max_r-padding, f"({cartesian_min_x:.2f}, {max_r:.2f})", color='r', alpha=0.5)
        
    for ii, f in enumerate(formulas):
        L = names[ii]
        if show_names:
            ax.text(hcList[ii] + .2*padding, yList[ii] - .2*padding, L, color='k',
                 fontsize=fontsize, rotation=rotation, alpha=0.9)
        if show_formula:
            L = f + ', \n' + names[ii]
            ax.text(hcList[ii] + .2*padding, yList[ii] - .2*padding, L, color='k',
                 fontsize=fontsize, rotation=rotation, alpha=0.9)
        
    for edge in pathway_edges: # cleaned edges to use
        ax.annotate("",
            xy=layout[edge[0]], xycoords='data',
            xytext=layout[edge[1]], textcoords='data',
            arrowprops=dict(arrowstyle="->", color="0.5", alpha=0.2,
                            shrinkA=5, shrinkB=5,
                            patchA=None, patchB=None,
                            connectionstyle="arc3,rad=-0.3",
                            ),
            )
    ax.set_rorigin(-1)
    for p in ax.spines.values():
        # inner, polar, start, end
        p.set_color('r')
        p.set_linewidth(0.1)
        
    median_y = np.median(yList)
    ax.set_rgrids([median_y])
    ax.grid(color='r', linewidth=0.05)
    # if xticks:
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xlabels])
    ax.text(1.5, max(yList) + .5, "H:C", ) # 
    if not yLabel_off:
        ax.text(1.97*np.pi, 1.6, ylabels[ydata], rotation=283)
    
    ax.set_title(title) 
    plt.savefig(outfile)


def assign_new_values(hcList, yList, cluttered, distance_min=0.05, iter=5,
                      xweight=0.02, rweight=0.01):
    '''
    cluttered : e.g. [18, 47, 49, 57, 69, 71, 73, 74, 93]
    add random moves to clustered, and replace values in hcList, yList
    distance_min : target min euclidean distance btw any two nodes
    '''
    new = sorted([(hcList[ii], yList[ii], ii) for ii in cluttered])
    
    for ii in range(iter):
        replacements = [new[0]]
        for x in new[1:]:
            if euclidean(x[:2], replacements[-1][:2]) > distance_min:
                replacements.append(x)
            else:
                replacements.append(
                    (x[0]+ xweight*np.random.random_sample(), 
                     x[1]+rweight*np.random.random_sample(), x[2])
                )
        new = sorted(replacements)
            
    for x in replacements:
        hcList[x[2]] = x[0]
        yList[x[2]] = x[1]
    return hcList, yList

def part_twins(hcList, yList, threshod=0.01, xspacer=0.01, rspacer=0.01):
    '''
    Separate nodes of identical coordinates
    '''
    def _match_(LL, y, threshod): 
        _r_ = False
        for x in LL:
            if euclidean(x, y) < threshod:
                _r_ = True
        return _r_
        
    LL = list(zip(hcList, yList))
    new = [LL[0]]
    for pair in LL[1:]:
        if _match_(new, pair, threshod):
            new.append((pair[0]+xspacer, pair[1]+rspacer))
        else:
            new.append(pair)
    return [x[0] for x in new], [x[1] for x in new]

def optimize_list_coordinates(hcList, yList, 
                              clustersize=2,
                              distance_cut=0.1, 
                              optimization_iter=5,
                              xspacer=0.01, rspacer=0.01):
    '''
    add spacer to cluttered nodes
    Example cluster:
        3.5241459360629888 0.26
        3.5241459360629888 0.26
        3.5241459360629888 0.26
        3.4270902040763933 0.2826086956521739
        3.4776400644860788 0.2708333333333333
        3.5241459360629888 0.26
        3.4270902040763933 0.2826086956521739
        3.4776400644860788 0.2708333333333333
        3.5241459360629888 0.26
    
    '''
    YM = pdist(list(zip(hcList, yList)), 'euclidean')
    ZM = linkage(YM, method='complete')
    metClus = fcluster(ZM, distance_cut, criterion='distance')
    # Compile clusters
    metClusDict = {}
    for ii in range(len(metClus)):
        if metClus[ii] in metClusDict:
            metClusDict[ metClus[ii] ].append(ii)
        else:
            metClusDict[ metClus[ii] ] = [ii]

    # print("number of clusters: ", len(metClusDict.keys()))
    for k,v in metClusDict.items():
        if len(v) >= clustersize:
            hcList, yList = assign_new_values(hcList, yList, v, 
                                              distance_min=0.05,
                                              iter=optimization_iter)
            
    hcList, yList = part_twins(hcList, yList, threshod=0.1, xspacer=xspacer, rspacer=rspacer)
    return hcList, yList
   

def populate_cpd_coordinates(cpdDict, formula='neutral_formula'):
    # full cpds
    allcpds = []
    for C in cpdDict.values():
        if C[formula]:
            try:
                f = C[formula]
                if 'C' in f and 'H' in f:
                    hc, oc, nc, nops_c, ccs_offset = get_kvli_ratios(f, 0) #
                    allcpds.append(
                        {
                            'id': C['id'],
                            'name': C['name'],
                            formula: f,
                            'ratio_H_C': hc,
                            'ratio_O_C': oc,
                            'ratio_N_C': nc,
                            'ratio_NOPS': nops_c,
                            'ccs_offset': ccs_offset,
                            'centrality': None, # future use
                        }
                    )
                else:
                    print(C, f)
            except KeyError:
                print(C)
                
    return allcpds


def plot_Cartesian_vkli_pathway(pathway_nodes, pathway_edges, 
                                ydata='ratio_NOPS',
                                show_names=True, show_formula=False, 
                                fontsize=5, rotation=30, padding=0.01,
                                width=8, height=8,   # inch
                                ylim=None, # or max value of yaxis
                                title='', outfile='cartesian_vk_plot.pdf', dpi=300
                                ):
    hcList = [x['ratio_H_C'] for x in pathway_nodes]
    yList = [x[ydata] for x in pathway_nodes]
    formulas = [x['neutral_formula'] for x in pathway_nodes]
    names = [x['name'] for x in pathway_nodes]
    ylabels = {'ratio_NOPS': 'NOPS:C', 'ratio_O_C': 'O:C', 'ratio_H_C': 'H:C'} # proper text for axis labels
    
    layout = dict((x['id'], (hcList[ii], yList[ii])) for ii,x in enumerate(pathway_nodes))
    
    _nodes = [x['id'] for x in pathway_nodes]
    #sizes = [len(fdict[x[0]]) * 3 for x in kv_list_gem]
    
    text_offset = padding
    plt.figure(figsize=(width, height))
    ax = plt.subplot()
    ax.scatter(hcList, yList, marker='o', facecolors='blue', 
               linewidths=.2, edgecolors='m', s=48, alpha=0.5,
                )
    for ii, f in enumerate(formulas):
        L = names[ii]
        if show_names:
            ax.text(hcList[ii] + text_offset/2, yList[ii] + text_offset, L, 
                 fontsize=fontsize, rotation=rotation, alpha=0.5)
        if show_formula:
            L = f + ', \n' + names[ii]
            ax.text(hcList[ii] + text_offset/2, yList[ii] + text_offset, L, 
                    fontsize=fontsize, rotation=rotation, alpha=0.5)
        
    for edge in pathway_edges:
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
        
    ax.set_xlabel("H:C ratio")
    if ylim:
        ax.set_ylim((0, ylim))
    ax.set_ylabel(ylabels[ydata])
    ax.set_title(title) 
    plt.savefig(outfile)


# experimenting 
def nx_fruchterman_reingold(pathway_edges, posDict,
    k=0.05, fixed=None, iterations=5, threshold=0.0001, dim=2, 
):
    '''
    Modified from networkx._fruchterman_reingold
    '''
    import networkx as nx
    
    G = nx.from_edgelist(pathway_edges)
    A = nx.to_scipy_sparse_array(G)
    nnodes = len(hcList)

    pos = np.zeros((len(G), dim)) 
    for i, n in enumerate(G):
        if n in posDict:
            pos[i] = np.asarray(posDict[n])

    # the initial "temperature"  is about .1 of domain area (=1x1)
    # this is the largest step allowed in the dynamics.
    # We need to calculate this in case our fixed positions force our domain
    # to be much bigger than 1x1
    t = 0.001
    # simple cooling scheme.
    # linearly step down by dt on each iteration so last iteration is size dt.
    dt = t / (iterations + 1)
    delta = np.zeros((pos.shape[0], pos.shape[0], pos.shape[1]), dtype=A.dtype)
    # the inscrutable (but fast) version
    # this is still O(V^2)
    # could use multilevel methods to speed this up significantly
    for iteration in range(iterations):
        
        print(pos[0])
        
        # matrix of difference between points
        delta = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
        # distance between points
        distance = np.linalg.norm(delta, axis=-1)
        # enforce minimum distance of 0.01
        np.clip(distance, 0.01, None, out=distance)
        # displacement "force"
        displacement = np.einsum(
            "ijk,ij->ik", delta, (k * k / distance**2 - A * distance / k)
        )
        # update positions
        length = np.linalg.norm(displacement, axis=-1)
        length = np.where(length < 0.01, 0.1, length)
        delta_pos = np.einsum("ij,i->ij", displacement, t / length)
        if fixed is not None:
            # don't change positions of fixed nodes
            delta_pos[fixed] = 0.0
        pos += delta_pos
        # cool temperature
        t -= dt
        if (np.linalg.norm(delta_pos) / nnodes) < threshold:
            break
    return pos

