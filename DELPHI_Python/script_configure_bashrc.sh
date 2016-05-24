#!/bin/bash

echo "Configuring bashrc to allow Python scripts to use the DELPHI library in this directory"
echo "(NOTE: WORKS ONLY IF YOU USE BASH SHELL!)"

echo $'\n# paths for DELPHI_Python' >> ~/.bashrc
echo "PYTHONPATH="`pwd`":\$PYTHONPATH" >> ~/.bashrc
echo "export PYTHONPATH" >> ~/.bashrc
echo "PATH="`pwd`":\$PATH" >> ~/.bashrc
echo "export PATH" >> ~/.bashrc

echo $'\nDone\n'

echo $'DON\'T FORGET TO DO:\nsource ~/.bashrc'

