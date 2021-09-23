#!/usr/bin/python3
"""Transfer files from ``jobs/prepared`` to ``jobs/run``.

By default, existing files are ignored,
and non-existing files are brought over.
Other options include "clobber".

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
import os
import os.path as osp
import shutil

#Modified from the python standard library original
def copytree(src, dst, symlinks=False):
  names = os.listdir(src)
  os.makedirs(dst)
  errors = []
  for name in names:
    srcname = os.path.join(src, name)
    dstname = os.path.join(dst, name)
    try:
      if symlinks and os.path.islink(srcname):
        linkto = os.readlink(srcname)
        os.symlink(linkto, dstname)
      elif os.path.isdir(srcname):
        copytree(srcname, dstname, symlinks)
      else:
        if not osp.isfile(dstname):
          shutil.copy2(srcname, dstname)
    except OSError as why:
      errors.append((srcname, dstname, str(why)))
    # catch the Error from the recursive copytree so that we can
    # continue with other files
    except Error as err:
      errors.extend(err.args[0])
  try:
    shutil.copystat(src, dst)
  except OSError as why:
    # can't copy file access times on Windows
    if why.winerror is None:
      errors.extend((src, dst, str(why)))
  if errors:
    raise Error(errors)


#Path to the "jobs" folder
jobs_dir=osp.abspath(osp.join(osp.split(__file__)[0],"../jobs"))
#Paths to the "prepared" and "run" folders
prep_dir=osp.join(jobs_dir,"prepared")
run_dir=osp.join(jobs_dir,"run")

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument("subpath",help="directory structure at the level below `prepared` and `run`")
parser.add_argument("--clobber",action="store_true",help="first erase the existing run folder, if any")
cmdline=parser.parse_args()

#Source and destination paths
srcdir=osp.join(prep_dir,cmdline.subpath)
assert osp.isdir(srcdir), "%s does not exist."%srcdir
destdir=osp.join(run_dir,cmdline.subpath)

#Erase existing destination, if requested
if cmdline.clobber and osp.isdir(destdir):
  shutil.rmtree(destdir)

#Do the copy
copytree(srcdir,destdir)