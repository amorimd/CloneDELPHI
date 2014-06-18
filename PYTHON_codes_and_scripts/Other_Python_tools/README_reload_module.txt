# how to reload a module under ipython (after it has been modified, for instance):
# example with code from 'Impedance.py'
reload(sys.modules['Impedance'])
from Impedance import *
