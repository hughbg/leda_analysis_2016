#!/usr/bin/env python

"""
plot_s2p.py

Script for plotting S2P files from anritsu VNA
"""

import os, sys
import numpy as np
import pylab as plt

def read_s2p(filename):
    """ Read an S2P file """
    
    ff = open(filename, 'r')
    dd = []
    
    # Read file
    lines = ff.readlines()
    
    # find last comment line (this is row headers)
    colnames_row = 0
    data_found = False
    for line in lines:
        if line.startswith('!') or line.startswith('#'):
            if not data_found:
                colnames_row += 1
            pass
        else:
            data_found = True
            dd.append(map(float, line.split()))
        
    colnames = np.array(lines[colnames_row-1].split())[1:]
    return colnames, np.array(dd)
    

def plot_s11_logmag(filename, c0='#333333', label=None, ls='solid'):
    """ Plot S11 log mag """
    colnames, data = read_s2p(filename)
    plt.plot(data[:, 0], data[:, 1], c=c0, label=label, ls=ls)
    plt.xlabel(colnames[0])
    plt.ylabel(colnames[1])
    plt.minorticks_on()
    
def plot_s21_logmag(filename, c0='#333333', label=None, atten=0):
    """ Plot S11 log mag """
    colnames, data = read_s2p(filename)
    plt.plot(data[:, 0], data[:, 3] + atten, c=c0, label=label)
    plt.xlabel(colnames[0])
    plt.ylabel(colnames[1])
    plt.minorticks_on()
    
def plot_s2p(filename, save_fig=True, show_fig=True, s11=True, c0='#cc0000', c1='#333333'):
    """ Plot S2P file """
    colnames, data = read_s2p(filename)
    
    if len(colnames) != 9:
        print "Unknown S2P format!"
    else:
        
        #plt.figure("S2P", figsize=(8,8))
        
        plt.suptitle(filename)
        if s11:
            plt.subplot(2,1,1)
        else:
            plt.subplot(2,2,1)
            
        plt.plot(data[:, 0], data[:, 1], c=c0)
        plt.xlabel(colnames[0])
        plt.ylabel(colnames[1])
        plt.minorticks_on()
        
        if s11:
            plt.subplot(2,1,2)
        else:
            plt.subplot(2,2,3)
        plt.plot(data[:, 0], data[:, 2], c=c1)
        plt.xlabel(colnames[0])
        plt.ylabel(colnames[2])
        plt.minorticks_on()
        
        if not s11:
            plt.subplot(2,2,2)
            plt.plot(data[:, 0], data[:, 3], c=c0)
            plt.xlabel(colnames[0])
            plt.ylabel(colnames[1])
            plt.minorticks_on()
            
            plt.subplot(2,2,4)
            plt.plot(data[:, 0], data[:, 4], c=c1)
            plt.xlabel(colnames[0])
            plt.ylabel(colnames[2])
            plt.minorticks_on()
        
        if save_fig:
            plt.savefig(filename + '.pdf')
        if show_fig:
            plt.show()   
        else:
            pass
            #plt.clf()     

def plot_logmag_multifile(s11_filename, s12_filename, s21_filename, s22_filename, 
                          save_fig=True, show_fig=True, title='', c0='#333333'):
    """ Plot S2P file """
    colnames, data = read_s2p(s11_filename)
    
    if len(colnames) != 9:
        print "Unknown S2P format!"
    else:
        
        #plt.figure("S2P", figsize=(8,8))
        
        plt.suptitle(title)
        plt.subplot(2,2,1)
        colnames, data = read_s2p(s11_filename)
        plt.plot(data[:, 0], data[:, 1], c=c0)
        plt.xlabel(colnames[0])
        plt.ylabel("S11 [dB]")
        plt.minorticks_on()
        
        plt.subplot(2,2,2)
        colnames, data = read_s2p(s12_filename)
        plt.plot(data[:, 0], data[:, 3], c=c0)
        plt.xlabel(colnames[0])
        plt.ylabel("S12 [dB]")
        plt.minorticks_on()

        plt.subplot(2,2,3)
        colnames, data = read_s2p(s21_filename)
        plt.plot(data[:, 0], data[:, 3] + 40, c=c0)
        plt.xlabel(colnames[0])
        plt.ylabel("S21 [dB]")
        plt.minorticks_on()
        
        plt.subplot(2,2,4)
        colnames, data = read_s2p(s22_filename)
        plt.plot(data[:, 0], data[:, 1], c=c0)
        plt.xlabel(colnames[0])
        plt.ylabel("S22 [dB]")
        plt.minorticks_on()
        
        if save_fig and title != '':
            plt.savefig(title + ".pdf")
        if show_fig:
            plt.show()   
        else:
            pass    

if __name__ == "__main__":
    try:
        fname = sys.argv[1]
        #fname = '2b2a-rj45-Ncal_allTERMoutput.s2p'
    except:
        print "USAGE: ./plot_s2p.py FILENAME"
        exit()
    
    plot_s2p(fname)