import sys
import os
import fnmatch
import subprocess

def crop_all_pdfs(mydir, outdir):
    """
    crop all pdf files in a folder and save with same name to another folder.
    """


    if not os.path.exists(mydir):
        os.makedirs(mydir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    all_pdf_files = [ f for f in os.listdir(mydir) if fnmatch.fnmatch(f, '*.pdf')]

    proc_set = set()    
    for pdf_file in all_pdf_files:
        proc_set.add(subprocess.Popen(['pdfcrop', os.path.join(mydir, pdf_file), os.path.join(outdir, pdf_file) ]) )

    # wait until all fmriresample finishes
    for p in proc_set:
        if p.poll() is None:
            p.wait()
            
if __name__ == '__main__':
    crop_all_pdfs(sys.argv[1], sys.argv[2])
