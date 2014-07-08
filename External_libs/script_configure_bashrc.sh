#!/bin/bash

echo "Configuring bashrc to allow Python scripts to use the python modules you installed"
echo "(NOTE: WORKS ONLY IF YOU USE BASH SHELL!)"
echo $'\n# Environment variables for local Python modules' >> ~/.bashrc
echo "export PYMOD='local'" >> ~/.bashrc

read -e -p $'Modify the path to the numpy "site-packages" folder if it is not this one:\n' -i "numpy-install/lib64/python2.6/site-packages" FILEPATH
if [[ ! -d `pwd`/$FILEPATH ]]; then
  echo "path not found !!"
  exit;
fi
echo "export PY_NUMPY="`pwd`/$FILEPATH >> ~/.bashrc

read -e -p $'Modify the path to the scipy "site-packages" folder if it is not this one:\n' -i "scipy-install/lib64/python2.6/site-packages" FILEPATH
if [[ ! -d `pwd`/$FILEPATH ]]; then
  echo "path not found !!"
  exit;
fi
echo "export PY_SCIPY="`pwd`/$FILEPATH >> ~/.bashrc

read -e -p $'Modify the path to the matplotlib "site-packages" folder if it is not this one:\n' -i "matplotlib-install/lib64/python2.6/site-packages" FILEPATH
if [[ ! -d `pwd`/$FILEPATH ]]; then
  echo "path not found !!"
  exit;
fi
echo "export PY_MATPL="`pwd`/$FILEPATH >> ~/.bashrc

echo $'DON\'T FORGET TO DO:\nsource ~/.bashrc'

