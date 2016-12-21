#!/usr/bin/env python
#
# clean SyEM example folder
# jcanton@mech.kth.se


import os, shutil

#------------------------------------------------------------------------------

# clean pipeMeshNek
try:
	os.remove('generateMesh')
	print('removing generateMesh')
except OSError:
	print('generateMesh not present')

try:
	shutil.rmtree('pipeMeshNek')
	print('removing pipeMeshNek')
except OSError:
	print('pipeMeshNek not present')

# clean Nek5000 svn
try:
	shutil.rmtree('nek5_svn')
	print('removing nek5_svn')
except OSError:
	print('nek5_svn not present')

try:
	os.rename('makenek.bak', 'makenek')
	print('restoring original makenek')
except OSError:
	print('makenek not present')

# clean simulation files
try:
	os.remove('base.rea')
	os.remove('base2d.rea')
	os.remove('pipe.map')
	os.remove('pipe.rea')
	os.remove('pipe.re2')
	os.remove('makefile')
	os.remove('avg.mod')
	os.remove('sem.mod')
	os.remove('pipe.f')
	os.remove('compiler.out')
	shutil.rmtree('obj')
	os.remove('nek5000')
	os.remove('SESSION.NAME')
	os.remove('pipe.sch')
	os.remove('sem_restart.txt')
	print('removing simulation files')
except OSError:
	print('simulation files not present')
