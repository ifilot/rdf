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
import progressbar

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
		print "Grabbing atoms for dataset %04i" % (len(snapshots)+1)
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

	print "Start calculating RDF"

	# collect sizes
	nr_bins = int(math.ceil(cutoff / binsize))
	nr_atoms = atoms.shape[0]

	# construct progress bar
	pbar = progressbar.ProgressBar(max_value=nr_atoms, redirect_stdout=True)

	# count distances
	distances = []
	for i in range(0, nr_atoms):
		pbar.update(i)
		for j in range(i+1, nr_atoms):
			d = numpy.linalg.norm(atoms[i] - atoms[j], 2)
			if d < cutoff:
				distances.append(d)
	pbar.update(nr_atoms)

	# construct rdf matrix
	rdf = numpy.zeros((3, nr_bins))

	# loop over bins
	for distance in distances:

		# do not count distances larger than cut-off
		if distance > cutoff:
			continue

		# calculate bin number
		bin_nr = int(math.floor(distance / binsize))

		rdf[1,bin_nr] += 1

	for i in range(0, nr_bins):
		r2 = i * binsize
		r1 = (i-1) * binsize
		rdf[2,i] = 4./3. * math.pi * (r2**3 - r1**3) # calculation shell volume

	# divide by volume
	for i in range(0, nr_bins):
		rdf[0,i] = (i + 0.5) * binsize
		rdf[1,i] /= rdf[2,i]

	print "Done calculating RDF"

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
filename = "%04i.png" % 1
print "Writing to %s" % filename
plt.savefig(filename)
