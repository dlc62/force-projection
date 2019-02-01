#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np

def get_contents(filename):
  with open(filename,'r') as f:
    fc = f.readlines()
  return fc

def write_contents(filename,data):
  with open(filename,'w') as f:
    f.writelines(data)

def get_section(contents,start_key,stop_key,start_offset=0,stop_offset=0):

  for i,line in enumerate(contents):
    if start_key in line:
      start_index = i + start_offset
      break
  for i in range(start_index,len(contents)):
    line = contents[i]
    if stop_key in line:
      stop_index = i+1 + stop_offset
      break
  data = contents[start_index:stop_index]
  return data

def get_geom_grad(data):

  geometry = get_section(data,'Standard Nuclear Orientation','Nuclear Repulsion Energy',-1,-1)
  gradient = get_section(data,'Cartesian Gradient','Attempting to Generate Delocalized',2,-2)
  raw_geometry = [list(map(float,line.split()[2:])) for line in geometry[4:-1]]
  raw_gradient = [list(map(float,line.split()[2:])) for line in gradient] 

  return np.asarray(raw_geometry),np.asarray(raw_gradient)
  
def project(forces,coords,atoms):

  projected_forces = []

  for atom1, atom2 in zip(atoms[:-1],atoms[1:]):
     
    bond = coords[atom1,:] - coords[atom2,:]
    force = forces[atom1,:] - forces[atom2,:]
    f_proj = bond.dot(force)/np.sqrt(bond.dot(bond))
    projected_forces.append(f_proj)
  
  return np.asarray(projected_forces)

def process_inputs(qcoutfile,ringfile):

  # get coords, force data
  base = qcoutfile.split(".")[0]
  data = get_contents(qcoutfile)
  coords,forces = get_geom_grad(data)

  # get list of atoms in ring (must be provided in bonded order sequence)
  atom_data = get_contents(ringfile)
  atoms_in_ring = []
  for line in atom_data:
    atoms_in_ring += list(map(int,line.split()))
  atom_indices = [i-1 for i in atoms_in_ring]
   
  # compute projected forces, and take sum of absolute values
  projected_forces = project(forces,coords,atom_indices)

  return np.sum(np.absolute(projected_forces))

#========================================================================

if __name__ == "__main__":

  if len(sys.argv) == 3:
    total_projected_force = process_inputs(sys.argv[1],sys.argv[2])
    print(total_projected_force)
  else:
    print("Usage: project_forces.py <Q-Chem gradient-containing output file> <atoms in ring data file>") 
    sys.exit()
