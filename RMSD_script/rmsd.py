#! /usr/bin/env python3
import os
import glob,re
from pymol import cmd

#local path to pdb files
directory = "/Users/neba/Desktop/CyclicaVsFDA/RMSD_script/*.pdb"

#Function to use pymol align API
def align_allFiles_to_allFiles(files=directory,cutoff=2,cycles=5,debug=1,full_matrix=0,method='align'):
  """
  Aligns all models in a list to all other models in the list

  usage:
    align_allfiles_to_allfiles [file_list][selection][cutoff=2][cycles=5][debug=0][full_matrix=0][method='align']

        where method can be align, super, cealign or rmscur

        where, cutoff and cycles are options passed to the align command.

    Setting debug=1 prints more information to the terminal or external GUI.
    Setting full_matrix=1 prints out the full symmetric matrix, rather than
    simply the top-half matrix

    Example:
      align_allFiles_to_allFiles files=name1 name2 name3 name4, full_matrix=1
  """
  file_list = glob.glob(files)
  file_list.sort()
  extension = re.compile( '(^.*[\/]|\.(pdb|ent|brk))' )

  cutoff = int(cutoff)
  cycles = int(cycles)
  full_matrix = int(full_matrix)
  debug=int(debug)

  object_list = []
  for filename in file_list:
    object_list.append(extension.sub('',filename))

  rmsd = {}
  rmsd_list = []
  for i in range(len(file_list)):
    obj_name1 = extension.sub('',file_list[i])
    cmd.load(file_list[i],obj_name1)
    for j in range(i+1,len(file_list)):
      obj_name2 = extension.sub('',file_list[j])
      cmd.load(file_list[j],obj_name2)

      rmsd_val = cmd.align('%s & %s' % (object_list[j],selection),
                    '%s & %s' % (object_list[i],selection),
                    cutoff=cutoff,cycles=cycles)

      rmsd.setdefault(obj_name1,{})[obj_name2] = rmsd_val[0]
      rmsd_list.append((rmsd_val[0],obj_name1,obj_name2))
      if debug:
        print("Alignment of %s to %s:" % (obj_name2,obj_name1))
        print("     Initial rmsd_val: %6.3f for %d atoms" % (rmsd_val[3],rmsd_val[4]))
        print("     Final rmsd_val: %6.3f for %d atoms after %d cycles\n" % (rmsd_val[0],rmsd_val[1],rmsd_val[2]))

      cmd.delete(obj_name2)
    cmd.delete(obj_name1)

  rmsd_list.sort()
# loop over dictionary and print out matrix of final rmsd_val values
  if debug:
    for object_name in object_list[:-1]:
      print("%s: %s" % (object_name,str(rmsd[object_name])))

    for r in rmsd_list:
      print("%6.3f  %s  %s" % r)

  print("%6s" % " ", end=' ')
  if full_matrix:
# fill in other half of matrix
    for i in range(len(object_list)):
      for j in range(i+1,len(object_list)):
        rmsd.setdefault(object_list[j],{})[object_list[i]] = rmsd[object_list[i]][object_list[j]]
      rmsd[object_list[i]][object_list[i]] = 0

    for i in range(len(rmsd)):
      print("%6s" % object_list[i], end=' ')
    print("")
    for i in range(len(object_list)):
      print("%6s" % object_list[i], end=' ')
      for j in range(len(object_list)):
        print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
      print("")
  else:
    for i in range(len(rmsd)):
      print("%6s" % object_list[i+1], end=' ')
    print("")
    for i in range(len(object_list)):
      print("%6s" % object_list[i], end=' ')
      for k in range(i):
        print("%6s" % " ", end=' ')
      for j in range(i+1,len(object_list)):
        print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
      print("")

cmd.extend('align_allFiles_to_allFiles',align_allFiles_to_allFiles)
