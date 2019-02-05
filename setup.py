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

# pull and compile
os.chdir('./pipeMeshNek')
subprocess.call('git pull', shell=True)
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
	subprocess.call('/bin/bash maketools genmap n2to3 reatore2', shell=True)
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


#------------------------------------------------------------------------------
# setup simulation
#
print('\n\tSetting up simulation...')

# generate mesh
try:
	assert os.path.isfile('base.rea')
except AssertionError:
	subprocess.call('./generateMesh', shell=True)

# convert mesh
try:
	assert os.path.isfile('pipe.re2')
except AssertionError:
	tfile = open('tmp.in', 'w')
	tfile.write('base\n')
	tfile.write('pipe\n')
	tfile.close()
	subprocess.call('./nek5_svn/bin/reatore2 < tmp.in', shell=True)
	os.remove('tmp.in')

# setup pipe.rea
if os.path.isfile('pipe.rea'):
	ofile = open('pipe.rea',     'r')
	nfile = open('pipe.rea.tmp', 'w')
	lines = ofile.readlines()
	for line in lines:
		if 'P002' in line:
			nfile.write('   -5300.0     P002: VISCOS\n')
		elif 'P011' in line:
			nfile.write('      1000     P011: NSTEPS\n')
		elif 'P015' in line:
			nfile.write('       100     P015: IOSTEP\n')
		elif 'P054' in line:
			nfile.write('   0.00000     P054: fixed flow rate dir: |p54|=1,2,3=x,y,z\n')
		elif 'P055' in line:
			nfile.write('   0.00000     P055: vol.flow rate (p54>0) or Ubar (p54<0)\n')
		elif 'P068' in line:
			nfile.write('         0     P068: iastep: freq for avg_all (0=iostep)\n')
		else:
			nfile.write(line)
	ofile.close()
	nfile.close()
	os.rename('pipe.rea.tmp', 'pipe.rea')

# generate map
try:
	assert os.path.isfile('pipe.map')
except AssertionError:
	tfile = open('tmp.in', 'w')
	tfile.write('pipe\n')
	tfile.write('0.05\n')
	tfile.close()
	subprocess.call('./nek5_svn/bin/genmap < tmp.in', shell=True)
	os.remove('tmp.in')

# compile nek
try:
	subprocess.call('./makenek pipe', shell=True)
except:
	print('\t*** Error compiling Nek5000 ***')

# write SESSION.NAME
try:
	assert os.path.isfile('SESSION.NAME')
except AssertionError:
	nfile = open('SESSION.NAME', 'w')
	nfile.write('pipe\n')
	nfile.write(cwd + '\n')
	nfile.close()

print('\tDone.\n')


#==============================================================================
print('\n---------')
print('All done.')
print('---------\n')
