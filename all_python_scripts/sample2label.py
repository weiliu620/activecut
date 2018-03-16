import os
import subprocess
import fnmatch
import sys
import socket

def sample2label(in_dir, out_dir):

    host_name = socket.gethostname()
    sample2label_set = set()
    
    for file in os.listdir(in_dir):
        sample2label_set.add( subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + host_name + '/sample2label', '-i', os.path.join(in_dir, file), '-o', os.path.join(out_dir, file)]) )

    for p in sample2label_set:
        if p.poll() is None:
            p.wait()



if __name__ == '__main__':
    sample2label(sys.argv[1], sys.argv[2])            
