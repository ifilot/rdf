#!/usr/bin/env python2

#
# Purpose: Calculate radial distribution function from a LAMMPS calculation
# 
# License:
# 
# calculate_rdf.py (c) by Ivo Filot <i.a.w.filot@tue.nl>
#
# calculate_rdf.py is licensed under a
# Creative Commons Attribution 4.0 International License.
#
# You should have received a copy of the license along with this
# work. If not, see <http://creativecommons.org/licenses/by/4.0/>.
# 

# import libraries
import numpy, math
import matplotlib.pyplot as plt
import re
import sys

cutoff = 7 # angstrom
filename = sys.argv[1]

def grab_atoms(fhandle, nr_atoms):
	"""
	@brief      grab atoms from file handler
	
	@param      fhandle   file reading handle (object of "open command")
	@param      nr_atoms  nr of atom lines
	
	@return     Nx3 matrix of atoms
	"""
	atoms = numpy.zeros([0,3])
	for i in range(0, nr_atoms):
		line = fhandle.readline()
		m = re.match(r"(\d+) ([\d.]+) ([\d.]+) ([\d.]+)", line)
		x = float(m.group(2));
		y = float(m.group(3));
		z = float(m.group(4));
		atoms = numpy.vstack([atoms, numpy.array([x,y,z])])
	return atoms

def read_atom_set_from_file(filename):
	"""
	@brief      Read for each snapshot the atoms
	
	@param      filename  The filename
	
	@return     array containing matrices of atoms
	"""
	f = open(filename, 'r')
	nr_atoms = int(f.readline())
	snapshots = []
	while nr_atoms > 0:
		f.readline() # skip line
		snapshots.append(grab_atoms(f, nr_atoms))
		line = f.readline()
		if line:
			nr_atoms = int(line)
		else:
			nr_atoms = 0
	return snapshots

def calculate_rdf(atoms, binsize, cutoff):
	"""
	@brief      Calculates the rdf.
	
	@param      atoms     The atoms
	@param      binsize   The binsize
	@param      cutoff    The cutoff
	
	@return     The rdf.
	"""

	# collect sizes
	nr_bins = int(math.ceil(cutoff / binsize))
	nr_atoms = atoms.shape[0]

	# count distances
	distances = []
	for i in range(0, nr_atoms):
		for j in range(i+1, nr_atoms):
			d = numpy.linalg.norm(atoms[i] - atoms[j], 2)
			if d < cutoff:
				distances.append(d)

	# sort the distances
	distances.sort()

	# construct rdf matrix
	rdf = numpy.zeros((3, nr_bins))

	# loop over bins
	idx = 0
	for i in range(0, nr_bins):
		rdf[0,i] = i * binsize
		r1 = (float(i)-0.5) * binsize
		r2 = (float(i) + 0.5) * binsize
		rdf[2,i] = 4./3. * math.pi * (r2**3 - r1**3) # calculation shell volume

		if r2 > distances[-1]:
			continue

		while(distances[idx] < r2):
			rdf[1,i] += 1
			idx += 1

	# divide by volume
	for i in range(0, nr_bins):
		rdf[1,i] /= rdf[2,i]

	return rdf

# collect all data from file
snapshots = read_atom_set_from_file(filename)
nr_snapshots = len(snapshots)
print "Read %i datasets" % nr_snapshots

# produce a rdf graph for each snapshot and store on the HD
rdf = calculate_rdf(snapshots[-1], 0.1, cutoff)
plt.plot(rdf[0], rdf[1], '-o')
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,200))
plt.savefig("%04i.png" % 0)
