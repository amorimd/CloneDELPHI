#!/bin/bash

echo "Configuring bashrc to allow Python scripts to use SUSSIX"
echo "(NOTE: WORKS ONLY IF YOU USE BASH SHELL!)"

echo $'\n# Python paths for SUSSIX' >> ~/.bashrc
echo "PYTHONPATH="`pwd`":"`pwd`"/Sussix:\$PYTHONPATH" >> ~/.bashrc
echo "export PYTHONPATH" >> ~/.bashrc

echo $'\nDone\n'

echo $'DON\'T FORGET TO DO:\nsource ~/.bashrc'

