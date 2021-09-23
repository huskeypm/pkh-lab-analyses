#!/usr/bin/python3
"""Get a summary of the status tags in the job folders.

by Tom Pace, 2020"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
import os
import os.path as osp
from collections import OrderedDict as odict
from datetime import datetime

#Constants
tagspec="__"
checklen=2*len(tagspec)
exclusions=['wham']

def istag(filename):
  "Return true if the file is a tag file"
  return (len(filename)>=checklen) and (filename[:len(tagspec)] == tagspec) and (filename[-len(tagspec):] == tagspec)

#Path to the "jobs" folder
jobs_dir=osp.abspath(osp.join(osp.split(__file__)[0],"../jobs"))
#Path to the "run" folder
run_dir=osp.join(jobs_dir,"run")

#Parse command line arguments
parser = argparse.ArgumentParser(description=globals()['__doc__'])
parser.add_argument("jobset",help="directory structure at the level below `run`")
parser.add_argument("outfile",nargs="?",default=None,help="optional output file name, defaults to timestamped file in set run folder")
parser.add_argument("--jobs_by_tag",action="store_true",help="to provide the listing of all jobs by tag")
cmdline=parser.parse_args()

#Get the job set folder
setdir=osp.join(run_dir,cmdline.jobset)
assert osp.isdir(setdir), "%s does not exist."%srcdir

#Get the output file path
if cmdline.outfile is None:
  ofsuffix=datetime.now().strftime("-%Y%m%d-%H%M%S")
  ofname="status{}.yaml".format(ofsuffix)
  ofpath=osp.join(setdir,ofname)
else:
  ofpath=osp.abspath(cmdline.outfile)

#Get the individual job folders
jobdirs=[itm for itm in os.scandir(setdir) if itm.is_dir() and itm.name not in exclusions]
jobdirs.sort(key=lambda d: d.name)

#Get status tags in each job folder
tags_by_job=odict()
for jd in jobdirs:
  tags_by_job[jd.name]=[]
  for root,dirs,files in os.walk(jd):
    tags_by_job[jd.name]+=[f for f in files if istag(f)]

#Get job folders for each tag
jobs_by_tag={}
for job,taglist in tags_by_job.items():
  for tag in taglist:
    if not tag in jobs_by_tag.keys():
      jobs_by_tag[tag]=[]
    jobs_by_tag[tag].append(job)

#Get count of jobs with each tag
counts_by_tag={}
for tag,joblist in jobs_by_tag.items():
  counts_by_tag[tag]=len(joblist)

#Write the output file
with open(ofpath,"w") as fp:
  #Total number of jobs
  fp.write("num_jobs: %d\n"%len(jobdirs))
  #Counts by tag
  fp.write("counts_by_tag:\n")
  for tag,count in counts_by_tag.items():
    fp.write("  %s: %d\n"%(tag,count))
  #Jobs by tag
  if cmdline.jobs_by_tag:
    fp.write("jobs_by_tag:\n")
    for tag,joblist in jobs_by_tag.items():
      fp.write("  %s:\n"%tag)
      for job in joblist:
        fp.write('    - "%s"\n'%job)
  #Tags by job
  fp.write("tags_by_job:\n")
  for job,taglist in tags_by_job.items():
    tagstr=", ".join(taglist)
    fp.write('  "%s": [%s]\n'%(job,tagstr))


