#!/usr/bin/env python
#
# clean SyEM example folder
# jcanton@mech.kth.se


import os, shutil

#------------------------------------------------------------------------------
#
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
