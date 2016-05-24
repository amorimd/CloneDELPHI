#!/usr/bin/python

# some basic physical constants

def electron_param():

    ''' gives basic electron parameters
    Output: e [C], mass [kg], speed of light [m/s], rest energy [J] '''
    
    e=1.602176487e-19; # elementary charge
    m0=9.10938e-31; # electron mass in kg
    c=299792458; # speed of light
    E0=m0*c**2 # rest energy

    return e,m0,c,E0
    
    
def proton_param():

    ''' gives basic proton parameters
    Output: e [C], mass [kg], speed of light [m/s], rest energy [J] '''
    
    e=1.602176487e-19; # elementary charge
    m0=1.6726216e-27; # proton mass in kg
    c=299792458; # speed of light
    E0=m0*c**2 # rest energy

    return e,m0,c,E0
    
