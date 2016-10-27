# -*- coding: utf-8 -*-
from __future__ import division
import os
import json
import math
import numpy
from pymatgen.electronic_structure.boltztrap import BoltztrapRunner, BoltztrapAnalyzer
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.electronic_structure.core import Spin
from pymatgen.matproj.rest import MPRester
from pymatgen.core.structure import Structure

#Output Settings
PRINT_DEBUG = True
PRINT_INFO = True
SAVE_INFO = True


with MPRester("xxxxxxxxxxxxxx") as mp:
	banddata = mp.get_data("mp-1186", prop="bandstructure")
	bs = banddata[0]["bandstructure"]
	structdata = mp.query("mp-1186", ["structure"])
	bs.structure = structdata[0]["structure"]

	i = 0
	while i < len(bs.kpoints):
		for band in bs.bands[Spin.up]:
			band.pop(i)
		bs.kpoints.pop(i)
		i = i + 1

	existkpoints = []
	i = 0
	while i < len(bs._kpoints):
		if any(numpy.array_equal(bs._kpoints[i].frac_coords, kpt) for kpt in existkpoints):
			for band in bs._bands[Spin.up]:
				band.pop(i)
			bs._kpoints.pop(i)
		else:
			existkpoints.append(bs._kpoints[i].frac_coords)
			i = i + 1
	nelec = 44 

	x = BoltztrapRunner(bs, nelec, dos_type="HISTO", energy_grid=0.0005, lpfac=5,
				run_type="BOLTZ", band_nb=None, tauref=0, tauexp=0, tauen=0,
				soc=False, doping=None, energy_span_around_fermi=0.5, scissor=0.0,
				kpt_line=None, spin=None, cond_band=False, tmax=1300.0, tgrid=50.0)
 	outfile = x.run(path_dir=os.getcwd(), convergence=True, write_input=True, clear_dir=False, max_lpfac=150, min_egrid=0.00005)
	y = BoltztrapAnalyzer.from_files(outfile)
#	EffMass = y.get_eig_average_eff_mass_tensor(temperature=300, doping=1e+18)

