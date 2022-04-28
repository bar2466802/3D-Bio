# Copyright (c) 2006, Sunjoong LEE <sunjoong@gmail.com>
# Licensed under the Academic Free License version 3.0
# http://sunjoong.sourceforge.net/
__author__ = 'Sunjoong LEE <sunjoong@gmail.com>'
__date__ = '2006-07-08'
__version__ = '1.0.0'
__all__ = ['TMalignWrap']

from itertools import izip
from getopt import gnu_getopt
from pyTMalign import pyTMalign, pyTMalignComb
from sys import exit


_pyTMalignLogo = """
 **************************************************************************
 *                                pyTMalign                               *
 * A python wrapper of Yang Zhang's TM-align program                      *
 * Comments on the warpper, please email to: sunjoong@gmail.com           *
 **************************************************************************"""

_pyTMalignResult = """
Chain 1:%-10s,  Size=%4d
Chain 2:%-10s,  Size=%4d (TM-score is normalized by %4d)

Aligned length=%4d, RMSD=%6.2f, TM-score=%7.5f, ID=%5.3f

(":" denotes the residue pairs of distance < 5.0 Angstrom)
%s
%s
%s
"""

_pyTMalignBriefHeader = """
Chain 1    Size Chain 2    Size Norm Aligned   RMSD TM-score    ID
------------------------------------------------------------------"""

_pyTMalignBriefResult = '%-10s %4d %-10s %4d %4d    %4d %6.2f  %7.5f %5.3f'

_pyTMalignUsage = """
 USAGE: %s [options] structure.pdb target.pdb

 OPTIONS:

        -h
        --help                Print this help message.
        -o        filename
        --outname=filename    Generate superimposed structure file.
        -a
        --ave                 Use average length.
        -s     sort_string
        --sort=sort_string    Sort by sort_string.
                              Default is no sorting.
                              sort_string example)
                                  l       Sort by aligned length
                                  r       Sort by rmsd
                                  s       Sort by TM-score
                                  i       Sort by sequence identity
                                  is      First sort key i and second s
                                  ir      First sort key i and second r
                                  il      First sort key i and second l
                                  isr     First sort key i and s and r
                                  isl     First sort key i and s and l
                                  isrl    First sort key i and s and r and l
                                  rsil    First sort key r and s and i and l
        -c
        --combination         Set combination mode.
        -b
        --brief               Set brief mode.
        -L      length
        --lengh=length        Use given length(integer value.)
        -d      cutoff
        -d0     cutoff
        --d0=cutoff           Set cutoff(real value.)
        -m      cutoff
        -dmin   cutoff
        --min=cutoff          Set minimum cutoff(real value.)

 EXAMPLES:

        $ %s structure.pdb target.pdb

        $ %s -o TM.sup structure.pdb target.pdb
        $ %s target.pdb -o TM.sup

        $ %s -a structure.pdb target.pdb
        $ %s structure.pdb target.pdb -a

        $ %s -a -o TM.sup structure.pdb target.pdb
        $ %s -a structure.pdb target.pdb -o TM.sup
        $ %s -o TM.sup structure.pdb target.pdb -a

        $ %s structure1.pdb structure2.pdb target.pdb
        $ %s -a -b -s s structure1.pdb structure2.pdb target.pdb

        $ %s -b pdb1.pdb pdb2.pdb pdb3.pdb pdb4.pdb
        $ %s -a -b -c -s i pdb1.pdb pdb2.pdb pdb3.pdb pdb4.pdb
"""

_TMalignLogo = """
 **************************************************************************
 *                                TM-align                                *
 * A protein structural alignment algorithm based on TM-score             *
 * Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9  *
 * Comments on the program, please email to: yzhang@ku.edu                *
 **************************************************************************"""

_TMalignUsage = """
 Brief instruction for running TM-align program:
 (For detail: Zhang & Skolnick, Nucl. Acid Res.2005 33, 2303)

 1. Align 'structure.pdb' to 'target.pdb'
   (By default, TM-score is normalized by the length of 'target.pdb')
   >TMalign structure.pdb target.pdb

 2. Run TM-align and output the superposition to 'TM.sup':
    (To view the superimposed structures by rasmol: >rasmol -script TM.sup)
   >TMalign structure.pdb target.pdb -o TM.sup

 3. If you want TM-score normalized by an assigned length, e.g. 100 aa:
   >TMalign structure.pdb target.pdb -L 100"""




def TMalignWrap(argv):
    """Wraper function of 'TMalign' command.

    This file, TMalignWrap.py could be executed in shell propmpt by
    $ python TMalignWrap.py structure.pdb target.pdb

    And then, this funcion will be called.

    See also:
    The pyTMalign() function's help.
    """
    structure, target, ave, sort, comb, brief, \
               outname, L_fix, d0_min, d0_fix = __TMalignWrap_getOpt(argv)

    if comb and type(structure) == type([]):
        aligned = pyTMalignComb(structure + [target],
                                ave, sort, outname, L_fix, d0_min, d0_fix)
    else:
        aligned = pyTMalign(structure, target,
                            ave, sort, outname, L_fix, d0_min, d0_fix)

    print _pyTMalignLogo
    print _TMalignLogo
    if type(structure) == type(''):
        structure_size, target_size, normalized_size, \
                        aligned_length, rmsd, tm_score, seq_id, \
                        aligned_seq1, aligned_seq2, aligned_seq3 = aligned
        if brief:
            print _pyTMalignBriefHeader
            print _pyTMalignBriefResult % (structure, structure_size,
                                           target, target_size,
                                           normalized_size, aligned_length,
                                           rmsd, tm_score, seq_id)
            print
        else:
            print _pyTMalignResult % (structure, structure_size,
                                      target, target_size, normalized_size,
                                      aligned_length, rmsd, tm_score, seq_id,
                                      aligned_seq1, aligned_seq3, aligned_seq2)
    elif comb:
        if brief:
            print _pyTMalignBriefHeader
            for (istructure, structure_size,
                 itarget, target_size, normalized_size,
                 aligned_length, rmsd, tm_score, seq_id,
                 aligned_seq1, aligned_seq2, aligned_seq3) in \
                 aligned.__iter__():
                print _pyTMalignBriefResult % (istructure, structure_size,
                                               itarget, target_size,
                                               normalized_size,
                                               aligned_length,
                                               rmsd, tm_score, seq_id)
            print
        else:
            for (istructure, structure_size,
                 itarget, target_size, normalized_size,
                 aligned_length, rmsd, tm_score, seq_id,
                 aligned_seq1, aligned_seq2, aligned_seq3) in \
                 aligned.__iter__():
                print _pyTMalignResult % (istructure, structure_size,
                                          itarget, target_size,
                                          normalized_size, aligned_length,
                                          rmsd, tm_score, seq_id,
                                          aligned_seq1, aligned_seq3,
                                          aligned_seq2)
    else:
        if brief:
            print _pyTMalignBriefHeader
            for (istructure, (structure_size, target_size, normalized_size,
                              aligned_length, rmsd, tm_score, seq_id,
                              aligned_seq1, aligned_seq2, aligned_seq3)) \
                              in izip(structure, aligned):
                print _pyTMalignBriefResult % (istructure, structure_size,
                                               target, target_size,
                                               normalized_size,
                                               aligned_length,
                                               rmsd, tm_score, seq_id)
            print
        else:
            for (istructure, (structure_size, target_size, normalized_size,
                              aligned_length, rmsd, tm_score, seq_id,
                              aligned_seq1, aligned_seq2, aligned_seq3)) \
                              in izip(structure, aligned):
                print _pyTMalignResult % (istructure, structure_size,
                                          target, target_size,
                                          normalized_size, aligned_length,
                                          rmsd, tm_score, seq_id,
                                          aligned_seq1, aligned_seq3,
                                          aligned_seq2)




def __TMalignWrap_getOpt(argv):
    name = argv[0]
    args = argv[1:]
    optlist, args = gnu_getopt(args, 'ho:L:m:d:as:cb',
                               ['help', 'outname=', 'length=', 'min=', 'd0=',
                                'ave', 'sort=', 'combination', 'brief'])
    optlist = dict(optlist)

    if '?' in args:
        index = args.index('?')
        optlist['-h'] = ''
        del args[index]
    if '-dmin' in args:
        index = args.index('-dmin')
        optlist['-m'] = args[index + 1]
        args[index:index + 2] = []
    if '-d0' in args:
        index = args.index('-d0')
        optlist['-d'] = args[index + 1]
        args[index:index + 2] = []

    opt = optlist.get('--help')
    if opt is not None:
        optlist['-h'] = opt
    opt = optlist.get('--outname')
    if opt:
        optlist['-o'] = opt
    opt = optlist.get('--length')
    if opt:
        optlist['-L'] = opt
    opt = optlist.get('--min')
    if opt:
        optlist['-m'] = opt
    opt = optlist.get('--d0')
    if opt:
        optlist['-d'] = opt
    opt = optlist.get('--ave')
    if opt is not None:
        optlist['-a'] = opt
    opt = optlist.get('--sort')
    if opt:
        optlist['-s'] = opt
    opt = optlist.get('--combination')
    if opt is not None:
        optlist['-c'] = opt
    opt = optlist.get('--brief')
    if opt is not None:
        optlist['-b'] = opt

    opt_error = False
    if optlist.has_key('-a'):
        ave = True
    else:
        ave = False
    if optlist.has_key('-s'):
        sort = optlist['-s']
        if len(sort) > 4:
            opt_error = True
        temp = list(sort)
        if 'l' in temp:
            del temp[temp.index('l')]
        if 'r' in temp:
            del temp[temp.index('r')]
        if 's' in temp:
            del temp[temp.index('s')]
        if 'i' in temp:
            del temp[temp.index('i')]
        if temp:
            opt_error = True
    else:
        sort = None
    if optlist.has_key('-c'):
        comb = True
    else:
        comb = False
    if optlist.has_key('-b'):
        brief = True
    else:
        brief = False
    if optlist.has_key('-o'):
        outname = optlist['-o']
    else:
        outname = None
    if optlist.has_key('-L'):
        if optlist['-L'].isalnum():
            L_fix = eval(optlist['-L'])
        else:
            opt_error = True
    else:
        L_fix = None
    if optlist.has_key('-m'):
        try:
            d0_min = float(optlist['-m'])
        except:
            opt_error = True
    else:
        d0_min = None
    if optlist.has_key('-d'):
        try:
            d0_fix = float(optlist['-d'])
        except:
            opt_error = True
    else:
        d0_fix = None

    if opt_error or optlist.has_key('-h') or len(args) < 2:
        print _TMalignUsage
        print _pyTMalignLogo
        print _pyTMalignUsage % (name, name, name, name, name,
                                 name, name, name, name, name,
                                 name, name, name)
        exit()

    if len(args) == 2:
        return (args[0], args[1],
                ave, sort, comb, brief, outname, L_fix, d0_min, d0_fix)
    else:
        return (args[:-1], args[-1],
                ave, sort, comb, brief, outname, L_fix, d0_min, d0_fix)




if __name__ == '__main__':
    from sys import argv
    TMalignWrap(argv)
