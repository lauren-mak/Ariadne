############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import with_statement
import os
import shlex
import shutil
import stat
import subprocess
import sys
import platform
import re
import gzip
import time
import random
from collections import defaultdict
from genericpath import isdir, exists

from os.path import join

import socket
socket.setdefaulttimeout(20.0)
from urllib2 import urlopen
import xml.etree.ElementTree as ET
import urllib

def download_fastq():
    dirs = os.listdir(main_dir)
    from joblib import Parallel, delayed
    Parallel(n_jobs=8)(delayed(download_one_fastq)(main_dir, dir)
             for dir in dirs)

def sra_id_to_sra_path(sra_id):
    ascp_pref = "" # /sra/sra-instant/reads/ByRun/sra/
    res = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/"+ sra_id[:3] +"/" + sra_id[:6] + "/"+ sra_id+ "/" + sra_id + ".sra"
    return res

def download_one_fastq(main_dir, sra_id):
    aspera_string = sra_id_to_sra_path(sra_id)
    if exists(join(main_dir, sra_id+"_1.fastq")):
        print sra_id + " already downloaded"
        return
    os.chdir(main_dir)
#    args = ['fastq-dump', '--split-files', sra_id]
    args = ['ascp', '-i', '/usr/local/etc/asperaweb_id_dsa.openssh', '-k', '1', '-T', '-l', '100m', aspera_string, '.']
#    print args
    print ('Downloading ' + sra_id)
    subprocess.call(args)
    print ('Downloaded ' + sra_id)
    return




def download_sra_from_list(sra_id_list):
    download_dir = sys.argv[2]
    sra_file = open(sra_id_list, 'r')
    for line in sra_file:
        download_one_fastq(download_dir, line.strip())


def main():
    if len(sys.argv) != 3:
        print "Usage: %s <list with SRA ids> <dir to download>" % sys.argv[0]
        sys.exit()
    sra_id_list = sys.argv[1]
    download_sra_from_list(sra_id_list)

if __name__ == '__main__':
    main()


