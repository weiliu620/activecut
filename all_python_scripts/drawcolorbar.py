from pylab import *
from numpy import *
import subprocess

def drawcolorbar(cmapname, out_file):
    """
    given colormap name, draw a colorbar and save it as a file
    """
    ioff()
    figure(figsize = (1, 10))
    axis("off")
    a = outer(arange(0, 1, 0.01), ones(10))
    imshow(a, aspect = 'auto', cmap = get_cmap(cmapname), origin = "lower")
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    savefig(out_file, dpi=100, bbox_inches='tight')
    subprocess.call(['/usr/bin/convert', '-trim', out_file, out_file ]) 
    
    


