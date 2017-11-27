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
from matplotlib import pyplot as plt

#
# assume Co here
#

fcc_unitcell = 2.50 * math.sqrt(2) * numpy.matrix([[0.5, 0.5, 0.0],[0.5, 0.0, 0.5],[0.0, 0.5, 0.5]])
fcc_atoms = numpy.matrix([0.0, 0.0, 0.0])
fcc = [fcc_unitcell, fcc_atoms, 'fcc']

bcc_unitcell = 2.50 * math.sqrt(2) * numpy.matrix([[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]])
bcc_atoms = numpy.matrix([0.0, 0.0, 0.0])
bcc = [bcc_unitcell, bcc_atoms, 'bcc']

hcp_unitcell = 2.50 * numpy.matrix([[1.0, 0.0, 0.0],[-0.5,math.sqrt(3.0)/2,0.0],[0.0,0.0,1.62]])
hcp_atoms = numpy.matrix([[0.0, 0.0, 0.0],[0.333,0.667,0.5]])
hcp = [hcp_unitcell, hcp_atoms, 'hcp']

sc_unitcell = 2.50 * numpy.matrix([[1.0, 0.0, 0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
sc_atoms = numpy.matrix([[0.0, 0.0, 0.0]])
sc = [sc_unitcell, sc_atoms, 'sc']

def calculate_rdf(struc, binsize, cutoff):
	"""
	@brief      Calculates the rdf.

	@param      unitcell  The unitcell
	@param      atoms     The atoms
	@param      binsize   The binsize
	@param      cutoff    The cutoff

	@return     The rdf.
	"""

	# expand
	[unitcell, atoms, title] = struc

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
		p1 = atoms[i].dot(unitcell)
		for j in range(i+1, nr_atoms):
			p2 = atoms[j].dot(unitcell)
			d = numpy.linalg.norm(p1 - p2, 2)
			distances.append(d)
		for j in range(0, nr_expanded_cell_atoms):
			p2 = expanded_cell_atoms[j].dot(unitcell)
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

	# divide by volume
	for i in range(0, nr_bins):
		rdf[1,i] /= rdf[2,i]

	return rdf


rdf_fcc = calculate_rdf(fcc, 0.01, 7)
rdf_hcp = calculate_rdf(hcp, 0.01, 7)
rdf_bcc = calculate_rdf(bcc, 0.01, 7)
rdf_sc = calculate_rdf(sc, 0.01, 7)

f, axarr = plt.subplots(2, 2)
axarr[0, 0].plot(rdf_fcc[0], rdf_fcc[1])
axarr[0, 0].set_title('FCC')
axarr[0, 0].set_xlabel('Distance r in Angstrom')
axarr[0, 0].set_ylabel('RDF g(r) [-]')
axarr[0, 1].plot(rdf_hcp[0], rdf_hcp[1])
axarr[0, 1].set_title('HCP')
axarr[0, 1].set_xlabel('Distance r in Angstrom')
axarr[0, 1].set_ylabel('RDF g(r) [-]')
axarr[1, 0].plot(rdf_bcc[0], rdf_bcc[1])
axarr[1, 0].set_title('BCC')
axarr[1, 0].set_xlabel('Distance r in Angstrom')
axarr[1, 0].set_ylabel('RDF g(r) [-]')
axarr[1, 1].plot(rdf_sc[0], rdf_sc[1])
axarr[1, 1].set_title('SC')
axarr[1, 1].set_xlabel('Distance r in Angstrom')
axarr[1, 1].set_ylabel('RDF g(r) [-]')
f.subplots_adjust(hspace=0.5)
plt.show()
