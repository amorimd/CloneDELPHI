#!/bin/bash

echo "Configuring bashrc to allow Python scripts to use DELPHI library"
echo "(NOTE: WORKS ONLY IF YOU USE BASH SHELL!)"

echo $'\n# library path for DELPHI' >> ~/.bashrc
echo "LD_LIBRARY_PATH="`pwd`":\$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH" >> ~/.bashrc

echo $'\nDone\n'

echo $'DON\'T FORGET TO DO:\nsource ~/.bashrc'

