import sys
import os
import fnmatch
import subprocess

def trim_all_pngs(indir, outdir):
    """
    crop all pdf files in a folder and save with same name to another folder.
    """


    if not os.path.exists(indir):
        os.makedirs(indir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    all_pdf_files = [ f for f in os.listdir(indir) if fnmatch.fnmatch(f, '*.png')]

    proc_set = set()    
    for pdf_file in all_pdf_files:
        proc_set.add(subprocess.Popen(['/usr/bin/convert', '-trim', os.path.join(indir, pdf_file), os.path.join(outdir, pdf_file) ]) )

    # wait until all fmriresample finishes
    for p in proc_set:
        if p.poll() is None:
            p.wait()

def crop_all_pngs(indir, outdir):
    """
    crop all pdf files in a folder and save with same name to another folder.
    """


    if not os.path.exists(indir):
        os.makedirs(indir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    all_pdf_files = [ f for f in os.listdir(indir) if fnmatch.fnmatch(f, '*.png')]

    proc_set = set()    
    for pdf_file in all_pdf_files:
        proc_set.add(subprocess.Popen(['/usr/bin/convert',os.path.join(indir, pdf_file),  '-crop', '800x315+0+140', os.path.join(outdir, pdf_file) ]) )

    # wait until all fmriresample finishes
    for p in proc_set:
        if p.poll() is None:
            p.wait()
                        
if __name__ == '__main__':
    crop_all_pngs(sys.argv[1], sys.argv[2])
