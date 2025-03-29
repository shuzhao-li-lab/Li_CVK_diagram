from .metModels import *
from .polarPlot import *

def run_demo(infile='../data/example_network_32.json'):
    net32 = json.load(open(infile))
    list_cpds, list_edges = net32['nodes'], net32['edges']
    cplot_LCVK_pathway(list_cpds, list_edges, 
                    formula='neutral_formula', cpd_name='name',
                      ydata='ratio_NOPS', 
                      min_theta=0.1*np.pi, max_theta=1.9*np.pi,
                      min_radius=0, max_radius=2.5,
                      padding=0.1, optimization=False,
                      optimization_spacer=0.3,
                      fontsize=5, rotation=30, show_names=True, show_formula=False, 
                      width=9, height=9,   # inch
                      title='Hello chemistry', outfile='demo_lcvk_plot.pdf', dpi=300
                      )

if __name__ == '__main__':
    run_demo()
