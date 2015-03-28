#!/usr/bin/env python

import tarfile
import urllib2
from os import makedirs
from os import path
from glob import glob
from BeautifulSoup import BeautifulSoup

# -- Some user defined parameters
ARCHIVE_DIR = 'Archives'

# -- Quantum Espresso archive webpage (2001-2009)
base_page_2001 = 'http://www.quantum-espresso.org/previous-versions/'
pre_path_2001 = ''

# -- Quantum Espresso archive webpages (2009-2012)
base_page_2009 = 'http://www.qe-forge.org/gf/project/q-e/frs/?action=' + \
                 'FrsReleaseBrowse&frs_package_id=18&_br_pkgrls_total' + \
                 '=50&_br_pkgrls_page=2'
pre_path_2009 = 'http://www.qe-forge.org'

# -- Quantum Espresso archive webpages (2012-2015)
base_page_2012 = 'http://www.qe-forge.org/gf/project/q-e/frs/?action=' + \
                 'FrsReleaseBrowse&frs_package_id=18'
pre_path_2012 = 'http://www.qe-forge.org'

# -- Gather all webpage addresses and pre paths
webpage_address_l = [base_page_2001, base_page_2009, base_page_2012]
pre_path_l = [pre_path_2001, pre_path_2009, pre_path_2012]


# -- Utility routines
def get_targz_webpage_links(
        webpage_address,
        pre_path='',
    ):
    response = urllib2.urlopen(webpage_address)
    html = response.read()
    soup = BeautifulSoup(html)
    targz_links = []
    for link in soup.findAll('a'):
        href = link.get('href')
        try:
            if href.endswith('.tar.gz') or href.endswith('.tgz'):
                targz_links += [pre_path + href]
        except:
            pass
    return targz_links

def get_all_targz_links(
        webpage_address_l,
        pre_path_l,
    ):
    all_links = []
    for address, pre_path in zip(webpage_address_l, pre_path_l):
        href_l = get_targz_webpage_links(address, pre_path)
        all_links += href_l
    return all_links

def filter_links(
        link_l,
        pattern_l=['espresso-', 'pw.'],
        exclude_pattern_l=['examples'],
    ):
    final_l = []
    for link in link_l:
        found_pattern = False
        for pattern in pattern_l:
            if pattern in link:
                found_pattern = True
                break
        for pattern in exclude_pattern_l:
            if pattern in link:
                found_pattern = False
                break
        if found_pattern:
            final_l += [link]
    return final_l

def download_archive(archive_address, base_dir=ARCHIVE_DIR):
    import urllib
    if not path.exists(base_dir):
        makedirs(base_dir)
    file_opener = urllib.URLopener()
    archive_name = path.basename(archive_address)
    file_archive_name = path.join(base_dir, archive_name)
    if not path.exists(file_archive_name):
        print 'Downloading: "%s"' % (archive_name,)
        try:
            file_opener.retrieve(archive_address, file_archive_name)
        except:
            print 'problem downloading archive...'
    return

# -- Untar and count number of lines in .f90 files
def get_total_line_number(archive_path):
    tar = tarfile.open(archive_path)
    all_fortran_files = [
                obj.name for obj in tar.getmembers()
                if obj.name.endswith('.f90')
            ]
    n_lines = 0
    for fortran_file in all_fortran_files:
        file_content = tar.extractfile(fortran_file).readlines()
        n_lines += len(file_content)
    return n_lines

# -- get all web addresses to PWscf or Espresso versions
all_links = get_all_targz_links(webpage_address_l, pre_path_l)
espresso_pwscf_links = filter_links(all_links)

# -- Download all archives (if there are not there already)
for link in espresso_pwscf_links:
    download_archive(link)

# -- Untar and count # of lines of Fortran90 code in archives
archive_fnames = sorted(glob(path.join(ARCHIVE_DIR, '*gz')))
print "%40s     %30s" % ('Archive name', '# of lines of F90',)
for archive_fname in archive_fnames:
    n_lines = get_total_line_number(archive_fname)
    print '%40s     %15i' % (archive_fname, n_lines,)
