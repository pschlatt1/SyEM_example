#!/usr/bin/env python
#
# setup SyEM example and fetch the necessary software
# jcanton@mech.kth.se


import os, shutil, subprocess, platform, stat

# get platform
pltfrm = platform.system()
# get current working directory
cwd = os.getcwd()


print('\nStarting setup for platform ' + pltfrm + '...\n')
#==============================================================================


#------------------------------------------------------------------------------
# import pipeMeshNek
#
print('\n\tImport pipeMeshNek...')

# clone
try:
	assert os.path.exists('./pipeMeshNek')
except AssertionError:
	subprocess.call('git clone https://github.com/jcanton/pipeMeshNek.git', shell=True)

# compile
try:
	assert os.path.isfile('generateMesh')
except AssertionError:
	os.chdir('./pipeMeshNek')
	subprocess.call('make', shell=True)
	os.rename('pipeMeshNek', '../generateMesh')
	os.chdir('..')

print('\tDone.\n')

#------------------------------------------------------------------------------
# import nek5000
#
print('\n\tImport Nek5000...')

# checkout revision 1093
try:
	assert os.path.exists('./nek5_svn')
except AssertionError:
	print("\n\t*** Accept nek's certificate ***\n")
	subprocess.call('svn co https://svn.mcs.anl.gov/repos/nek5 nek5_svn', shell=True)
	os.chdir('./nek5_svn')
	subprocess.call('svn up -r 1093', shell=True)
	os.chdir('..')

# compile tools
try:
	assert os.path.isfile('nek5_svn/bin/genmap')
except AssertionError:
	os.chdir('./nek5_svn/trunk/tools')
	ofile = open('maketools',     'r')
	nfile = open('maketools.tmp', 'w')
	lines = ofile.readlines()
	for line in lines:
		if 'mkdir -p' in line:
			nfile.write('mkdir -p ' + cwd + '/nek5_svn/bin\n')
		elif 'HOME' in line:
			nfile.write('bin_nek_tools="' + cwd + '/nek5_svn/bin"\n')
		elif 'mcmodel' in line:
			if (pltfrm == 'Darwin'):
				nline = line.replace('-mcmodel=medium', '')
				nfile.write(nline)
		else:
			nfile.write(line)
	ofile.close()
	nfile.close()
	os.rename('maketools.tmp', 'maketools')
	subprocess.call('/bin/bash maketools all', shell=True)
	os.chdir('../../../')

# we do not need this as we use a modified makenek
## # copy makenek
## try:
## 	assert os.path.isfile('makenek')
## except AssertionError:
## 	shutil.copy('nek5_svn/trunk/nek/makenek', 'makenek')
if os.path.isfile('makenek'):
	ofile = open('makenek',     'r')
	nfile = open('makenek.tmp', 'w')
	lines = ofile.readlines()
	for line in lines:
		if '$HOME' in line:
			nline = line.replace('$HOME', cwd)
			nfile.write(nline)
		else:
			nfile.write(line)
	ofile.close()
	nfile.close()
	os.rename('makenek',     'makenek.bak')
	os.rename('makenek.tmp', 'makenek')
	st = os.stat('makenek')
	os.chmod('makenek', st.st_mode | stat.S_IEXEC)


print('\tDone.\n')

#==============================================================================
print('Setup done.')
