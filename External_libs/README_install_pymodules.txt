#############################################
# INSTALLATION of Python modules locally
# (for instance on lxplus at CERN)
#############################################


#########
# NUMPY
#########

tar -xzvf numpy-1.7.0.tar.gz
cd numpy-1.7.0

# if necessary clean out
rm -rf build

# then (to install it in External_libs/numpy-install)

python setup.py build --fcompiler=gnu95
python setup.py install --prefix ../numpy-install

#########
# SCIPY
#########

tar -xzvf scipy-0.12.0.tar.gz
cd scipy-0.12.0

# change setup.py: add 2 lines (at line 145, in function setup_package) (to point to correct numpy installation directory)
import sys
sys.path.insert(1,'../numpy-install/lib64/python2.6/site-packages');
# NOTE: in the path above you might need to change "python2.6" with your 
# actual python version

# if necessary clean out
rm -rf build

# then (to install it in External_libs/scipy-install)
python setup.py build --fcompiler=gnu95
python setup.py install --prefix=../scipy-install

##############
# MATPLOTLIB
##############

tar -xzvf matplotlib-1.2.1.tar.gz
cd matplotlib-1.2.1

# change setupext.py: add 1 line (at line 67) (to point to correct numpy installation directory)
sys.path.insert(1,'../numpy-install/lib64/python2.6/site-packages');
# NOTE: in the path above you might need to change "python2.6" with your 
# actual python version

# if necessary clean out
rm -rf build

# create setup.cfg file (to avoid using Tkinter)
cp setup.cfg.template setup.cfg
# modify setup.cfg: uncomment the following lines (to avoid using Tkinter)
gtk = False
gtkagg = False
tkagg = False
macosx = False


# then (to install it in External_libs/matplotlib-install)
# (not sure the first make is useful)
make
python setup.py build
python setup.py install --prefix=../matplotlib-install

