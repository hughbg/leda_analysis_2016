#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import subprocess

def change_nec_freq(necname, freq):
    """Modify the FR card in a NEC input file to run at freq."""
    
    f = open(necname,'r+')
    buf = f.read()
    lines = buf.splitlines(True)
    # Substitute the freq in the right field of the FR card
    for i in range(len(lines)):
        if lines[i][:2] == 'FR':
            vals = re.split(',| +',lines[i])
            vals[5] = "%.2f" % freq
            lines[i] = " ".join(vals)
            # Make sure this line ends in newline
            if lines[i][-1] != '\n':
                lines[i] += '\n'
    # Rewrite the file
    f.seek(0)
    f.writelines(lines)
    f.truncate()
    f.close()


infiles = ['leda_xep_1.nec', 'leda_xet_1.nec', 'leda_yep_1.nec', 'leda_yet_1.nec']
for freq in [5, 10, 12, 15, 20, 25, 30, 35, 38, 40, 45, 50, 55, 60, 65, 70, 74, 75, 80, 85, 88, 90, 95, 100]:
    t = []
    for infile in infiles:
        outfile = infile.replace('_1', '_%i' % freq)
        os.system('cp %s %s' % (infile, outfile))
        change_nec_freq(outfile, freq)
        
        necfile = outfile.replace('.nec', '.out')
        t.append( subprocess.Popen(['nec4d', outfile, necfile]) )
    for p in t:
        p.wait()
        
