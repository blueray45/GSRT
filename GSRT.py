# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 03:03:06 2015

@author: irnakat
"""

import click
import IOfile
from TFCalculator import TFCalculator as TFC
import pylab as plt
import numpy as np

@click.command()
@click.option('--filename', default='input.dat', help='Name of configuration file.')
@click.option('--mode', default='sh-linear-elastic', \
    help='mode : sh-linear-elastic,sh-linear-viscoelastic')
@click.option('--method', default='tf_kramer286_sh', \
    help='mode : tf_kramer286_sh,tf_knopoff_sh,tf_knopoff_sh_adv')

def run_analysis(filename,mode,method):
    click.echo('Reading file : %s'%filename)
    data = IOfile.parsing_input_file(filename)
    click.echo('Creating class...')
    theclass = TFC(data)
    click.echo('Calculating transfer function using %s method'%method)
    if method=='tf_kramer286_sh':
        theclass.tf_kramer286_sh()
    elif method=='tf_knopoff_sh':
        theclass.tf_knopoff_sh()
    elif method=='tf_knopoff_sh_adv':
        theclass.tf_knopoff_sh_adv()
        
    plt.plot(theclass.freq,np.abs(theclass.tf[0]),label=method)
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Amplification')
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True,which='both')
    plt.legend(loc='best',fancybox=True,framealpha=0.5)
    #plt.axis('tight')
    plt.autoscale(True,axis='x',tight=True)
    plt.tight_layout()
    plt.savefig('test.png', format='png')
    click.echo(click.style('Calculation has been finished!',fg='green'))

if __name__ == '__main__':
    run_analysis()