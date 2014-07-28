#!/bin/bash

echo "Configuring bashrc to add this directory in your PATH"
echo "(NOTE: WORKS ONLY IF YOU USE BASH SHELL!)"

echo $'\n# paths for Data_analysis' >> ~/.bashrc
echo "PATH="`pwd`":\$PATH" >> ~/.bashrc
echo "export PATH" >> ~/.bashrc

echo $'\nDone\n'

echo $'DON\'T FORGET TO DO:\nsource ~/.bashrc'

