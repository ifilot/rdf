#!/usr/bin/env python2

#
# Purpose: Calculate theoretical Radial Distribution Function (RDF) for an FCC unit cell
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

import numpy, math
from matplotlib.pyplot import plot, show

# Ni bulk FCC
lattice_constant = 3.52

fcc_unitcell = lattice_constant * numpy.matrix([[0.5, 0.5, 0.0],[0.5, 0.0, 0.5],[0.0, 0.5, 0.5]])
fcc_atoms = numpy.matrix([0.0, 0.0, 0.0])

bcc_unitcell = lattice_constant * numpy.matrix([[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]])
bcc_atoms = numpy.matrix([0.0, 0.0, 0.0])

cutoff = 10 # angstrom

def calculate_rdf(unitcell, atoms, binsize, cutoff):
	"""
	@brief      Calculates the rdf.
	
	@param      unitcell  The unitcell
	@param      atoms     The atoms
	@param      binsize   The binsize
	@param      cutoff    The cutoff
	
	@return     The rdf.
	"""

	# collect sizes
	nr_bins = int(math.ceil(cutoff / binsize))
	nr_atoms = atoms.shape[0]

	# expand unit cell
	expanded_cell_atoms = numpy.zeros([0,3])
	dp = int(math.ceil(cutoff / numpy.linalg.norm(unitcell[0])))
	for p in range(-dp, dp):
		dq = int(math.ceil(cutoff / numpy.linalg.norm(unitcell[1])))
		for q in range(-dq, dq):
			dr = int(math.ceil(cutoff / numpy.linalg.norm(unitcell[2])))
			for r in range(-dr, dr):
				if p is 0 and q is 0 and r is 0:
					continue
				for i in range(0, nr_atoms):
					expanded_cell_atoms = numpy.vstack([expanded_cell_atoms, atoms[i]+numpy.array([p,q,r])])

	# count distances
	nr_expanded_cell_atoms = expanded_cell_atoms.shape[0]
	distances = []
	for i in range(0, nr_atoms):
		p1 = unitcell.dot(atoms[i].transpose())
		for j in range(0, nr_expanded_cell_atoms):
			p2 = unitcell.dot(expanded_cell_atoms[j].transpose())
			d = numpy.linalg.norm(p1 - p2, 2)
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

		while(distances[idx] < r2):
			rdf[1,i] += 1
			idx += 1

	# normalize for r0
	for i in range(0, nr_bins):
		rdf[1,i] /= rdf[2,i]

	# find value for first peak (for g0 normalization)
	idx = 0
	g0 = r0 = 0.0
	while rdf[1,idx] < 0.1:
		idx += 1
		g0 = rdf[1,idx]
		r0 = rdf[0,idx]

	# normalize for g0
	for i in range(0, nr_bins):
		rdf[1,i] /= g0
		rdf[0,i] /= r0

	plot(rdf[0], rdf[1], '-o')
	show()

calculate_rdf(fcc_unitcell, fcc_atoms, 0.1, cutoff)