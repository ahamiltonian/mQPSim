This repository is a Geant4 simulation of the MoEDAL mini-charged particle detector (mQP).

I've documented here what I did to install Geant4 and run the mQPSim code on my macbook running OSX 10.13.  


GEANT4 Installation
-------------------

1. ran 'xcode-select --install' to get clang installed

2. downloaded CMake from cmake.org and installed it
  - made it accessible to command line by running
    'sudo "/Applications/CMake.app/Contents/bin/cmake-gui" --install'

3. Installed XQuartz for macOS from https://www.xquartz.org/ (this gives me X11)

4. Download the data from http://geant4.web.cern.ch/support/download


(I tried the pre-compiled version of Geant4, but had some problems, so moved to installing from source code.  A Stack Overflow post that helped me:
https://stackoverflow.com/questions/23917587/geant4-example-b1-build-error-no-rule-to-make-target)

5. download source code from http://geant4.web.cern.ch/support/download
  - you should get a tarball named geant4.10.04.p01.tar.gz
  - unpack the tarball (just double click on the file)
  - you should get a folder named geant4.10.04.p01
  - copy that folder to somewhere that you want to have geant installed, I just moved it to my home directory /Users/andrew

6. Build Geant4
  - make a directory called build (ie. do 'mkdir build')
  - move to that directory (ie. do 'cd build')
  - do 'cmake -DCMAKE_INSTALL_PREFIX=../geant4-install -DGEANT4_USE_QT=ON ../geant4.10.04.p01'
  - copy all the data files from http://geant4.web.cern.ch/support/download into /Users/andrew/geant4-install/share/Geant4-10.4.1/data/
  - do 'make -j2'  (change 2 to 4 if you have a quad core processor, it will be faster)  this step can take an hour or so, and needs to download about 2GB of data
  - do 'make install'
  - do 'cd ../geant4-install/bin/'
  - do 'source ./geant4.sh'

7. Run a Geant4 example to be sure Geant4 is working
  - go to the examples directory 'cd ~/geant4-install/share/Geant4-10.4.1/examples/basic/B1/'
  - make a build directory, then go there (ie. do 'mkdir build' then 'cd build')
  - do 'cmake ../'
  - do 'make'

8. Download and run mQPSim
  - start a github account https://github.com/join
  - download the mQPSim code by running 'git clone https://github.com/ahamiltonian/mQPSim.git'
  - this will create an mQPSim directory, go into that directory
  - make a build directory and move into it
  - run 'cmake ../'
  - run 'make', this should create an executable named mQPSim
  - test the executable by running './mQPSim'

** Notes on the graphics viewer
- they all seem outdated...
- Wt looks like it is actually up to date, I'm trying that one...
- https://www.webtoolkit.eu/wt/documentation
   - needed to get boost libraries https://www.boost.org/users/download/
   - frustrated...  could not get wt to see the boost libraries
   - aborting...


** When you open a new terminal, you may need to do 'source ~/geant4-install/bin/geant4.sh'


ROOT Installation
------------------

1. go to root.cern.ch -> download -> latest "pro" release ("pro" is for production, not professional)
2. I downloaded the dmg for OSX 10.13
3. run the dmg to install
4. add root to your command line by adding "source /Applications/root_v6.12.06/bin/thisroot.sh" to your .bash_profile


WHAT TO DO NEXT LIST...
-----------------------

- added hits collection
- need to make light bounce in scintillator
  - tryed to use the Goddess package (https://forge.physik.rwth-aachen.de/projects/goddess-package/wiki), could not get boost library to work properly, goddess could not find boost header files, frustrated...
  -
