#!/usr/bin/env python

try:
    from numpy.distutils.core import setup, Extension
except:
    from sys import exit
    print 'Error: Install \'NumPy\' and \'F2PY\' first!'
    exit()
try:
    from sunjoong.FileUtils import openByFileExtension
except:
    from sys import exit
    print 'Error: Install \'Sunjoong\' first!'
    exit()


name = 'pyTMalign'
version = '1.0'
author = 'Sunjoong LEE'
author_email = 'sunjoong@gmail.com'
description = 'Python wrapper of Yang Zhang\'s TM-align and more'
long_description = """
Python wrapper of Yang Zhang's TM-align and more

Example1)
"pyTMalign" script is equivalent to "TMalign" program;

$ pyTMalign 153l.pdb 1dps.pdb

Example2)
"pyTMalign" script can read ".gz", ".bz2" and ".Z" files;

$ pyTMalign 153l.pdb.gz 1dps.pdb.Z

Example3)
"pyTMalign" script supports position independent options;

$ pyTMalign -a 153l.pdb 1dps.pdb
$ pyTMalign 153l.pdb 1dps.pdb -a

Example4)
"pyTMalign" script can take more than 2 files to align;

$ pyTMalign 153l.pdb 1e70.pdb 1qp8.pdb 1dps.pdb

"TMalign" program takes only 2 files and above example is like that;

$ TMalign 153l.pdb 1dps.pdb
$ TMalign 1e70.pdb 1dps.pdb
$ TMalign 1qp8.pdb 1dps.pdb

Example5)
The sort option could be set in "pyTMalign" script;

$ pyTMalign -s s 153l.pdb 1e70.pdb 1qp8.pdb 1dps.pdb

Example6)
The combination option is supported;

$ pyTMalign -b 153l.pdb 1e70.pdb 1qp8.pdb 1dps.pdb

It is like that in "TMalign" program;

$ TMalign 153l.pdb 1dps.pdb
$ TMalign 1e70.pdb 1dps.pdb
$ TMalign 1qp8.pdb 1dps.pdb
$ TMalign 153l.pdb 1qp8.pdb
$ TMalign 1e70.pdb 1qp8.pdb
$ TMalgin 153l.pdb 1e70.pdb

Example7)
The brief option is supported;

$ pyTMalign -a -b -c -s i 153l.pdb 1dps.pdb 1e70.pdb 1qp8.pdb 1xdw.pdb 2fcf.pdb

Bellow is it's output;

 **************************************************************************
 *                                pyTMalign                               *
 * A python wrapper of Yang Zhang's TM-align program                      *
 * Comments on the warpper, please email to: sunjoong@gmail.com           *
 **************************************************************************

 **************************************************************************
 *                                TM-align                                *
 * A protein structural alignment algorithm based on TM-score             *
 * Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9  *
 * Comments on the program, please email to: yzhang@ku.edu                *
 **************************************************************************

Chain 1    Size Chain 2    Size Norm Aligned   RMSD TM-score    ID
------------------------------------------------------------------
1qp8.pdb    295 2fcf.pdb     89  192      52   4.14  0.19048 0.016
153l.pdb    185 1e70.pdb    499  342      81   5.82  0.15794 0.040
1dps.pdb    159 1xdw.pdb    331  245      79   5.21  0.20907 0.042
153l.pdb    185 1dps.pdb    159  172      69   5.40  0.25623 0.042
153l.pdb    185 1qp8.pdb    295  240      96   6.04  0.22797 0.044
1e70.pdb    499 1qp8.pdb    295  397     149   6.13  0.24539 0.047
1e70.pdb    499 1xdw.pdb    331  415     175   6.58  0.26443 0.047
1dps.pdb    159 1e70.pdb    499  329      99   4.77  0.22097 0.049
1dps.pdb    159 2fcf.pdb     89  124      38   4.65  0.17406 0.049
1dps.pdb    159 1qp8.pdb    295  227      84   5.39  0.22999 0.050
153l.pdb    185 2fcf.pdb     89  137      45   4.22  0.21298 0.052
153l.pdb    185 1xdw.pdb    331  258      92   5.88  0.21305 0.067
1e70.pdb    499 2fcf.pdb     89  294      54   4.80  0.13005 0.070
1xdw.pdb    331 2fcf.pdb     89  210      56   4.69  0.17531 0.091
1qp8.pdb    295 1xdw.pdb    331  313     283   3.52  0.73923 0.208

Example8)
You can use python wrapper module within your python code.
You might be interested in pyTMalign() function.

>>> from sunjoong.pyTMalign import pyTMalign
>>>
>>> structure = '153l.pdb'
>>> target = '1dps.pdb'
>>>
>>> structure_size, target_size, normalized_size, aligned_length, \\
... rmsd, tm_score, seq_id, \\
... aligned_seq1, aligned_seq2, aligned_seq3 = \\
... pyTMalign(structure, target)
>>>
>>> print structure_size, target_size, normalized_size, aligned_length
>>> print rmsd, tm_score, seq_id
>>> print aligned_seq1
>>> print aligned_seq2
>>> print aligned_seq3

Example9)

>>> size1, size2, normalized_size, aligned_length, \\
... rmsd, tm_score, seq_id, \\
... aligned_seq1, aligned_seq2, aligned_seq3 = \\
... pyTMalign(structure, target, ave = True)
>>>

Example10)

>>> structure1 = '153l.pdb'
>>> structure2 = '1e70.pdb'
>>> structure3 = '1qp8.pdb'
>>> structure = [structure1, structure2, structure3]
>>>
>>> target = '1dps.pdb'
>>>
>>> aligned = pyTMalign(structure, target)
>>>
>>> structure_size, target_size, normalized_size, aligned_length, \\
... rmsd, tm_score, seq_id, \\
... aligned_seq1, aligned_seq2, aligned_seq3 = aligned[0]
>>> print structure[0], structure_size
>>> print target, target_size, normalized_size
>>> print normalized_size, aligned_length, rmsd,tm_score, seq_id
>>>
>>> structure_size, target_size, normalized_size, aligned_length, \\
... rmsd, tm_score, seq_id, \\
... aligned_seq1, aligned_seq2, aligned_seq3 = aligned[1]
>>> print structure[1], structure_size
>>> print target, target_size, normalized_size
>>> print normalized_size, aligned_length, rmsd,tm_score, seq_id
>>>

Example11)

>>> aligned = pyTMalign(structure, target, ave = True, sort = 'i')

Example12)

>>> from sunjoong.pyTMalign import pyTMalignComb
>>>
>>> structure1 = '153l.pdb'
>>> structure2 = '1e70.pdb'
>>> structure3 = '1qp8.pdb'
>>> structure4 = '1dps.pdb'
>>> structure = [structure1, structure2, structure3, structure4]
>>>
>>> aligned = pyTMalignComb(structure)
>>>
>>> istructure, structure_size, target, target_size, normalized_size, \\
... aligned_length, rmsd, tm_score, seq_id, \\
... aligned_seq1, aligned_seq2, aligned_seq3 = aligned[0]
>>> print istructure, structure_size
>>> print itarget, target_size, normalized_size
>>> print aligned_length, rmsd, tm_score, seq_id
>>>
>>> istructure, structure_size, target, target_size, normalized_size, \\
... aligned_length, rmsd, tm_score, seq_id, \\
... aligned_seq1, aligned_seq2, aligned_seq3 = aligned[1]
>>> print istructure, structure_size
>>> print itarget, target_size, normalized_size
>>> print aligned_length, rmsd, tm_score, seq_id
>>>

Example13)

>>> aligned = pyTMalignComb(structure, ave = True, sort = 'i')

Example14)

>>> from sunjoong.pyTMalign import pyTMalignMatrix
>>>
>>> structure1 = '153l.pdb'
>>> structure2 = '1e70.pdb'
>>> structure3 = '1qp8.pdb'
>>> structure4 = '1dps.pdb'
>>> structure = [structure1, structure2, structure3, structure4]
>>>
>>> aligned = pyTMalignMatrix(structure)
>>>
>>> istructure, structure_size, target, target_size, normalized_size, \\
... aligned_length, rmsd, tm_score, seq_id, \\
... aligned_seq1, aligned_seq2, aligned_seq3 = aligned.next()
>>> print istructure, structure_size
>>> print itarget, target_size, normalized_size
>>> print aligned_length, rmsd, tm_score, seq_id
>>>
>>> istructure, structure_size, target, target_size, normalized_size, \\
... aligned_length, rmsd, tm_score, seq_id, \\
... aligned_seq1, aligned_seq2, aligned_seq3 = aligned.next()
>>> print istructure, structure_size
>>> print itarget, target_size, normalized_size
>>> print aligned_length, rmsd, tm_score, seq_id
>>>


Special thanks to my mentor, Prof. Dr. Jooyoung Lee.
I saw he used "TMalign" at his work.
So, I knew there is a good tool of "TMalign" and wanted to improve it.

I thanks to Prof. Dr. Yang Zhang, too.
He is the author of the orignal "TMalign" program
and gave me quick response of my questions.


2006.07.08    Sunjoong LEE <sunjoong@gmail.com>

"""
url = 'http://sunjoong.sourceforge.net/'
license = 'Academic Free License 3.0'
classifiers = filter(None, """
License :: OSI-Approved Open Source :: Academic Free License (AFL)
Intended Audience :: by Industry or Sector :: Science/Research
Development Status :: 5 - Production/Stable
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Fortran
Programming Language :: Python
Operating System :: Grouping and Descriptive Categories :: All POSIX (Linux/BSD/UNIX-like OSes)
User Interface :: Textual :: Command-line
""".split('\n'))


py_modules = ['sunjoong.pyTMalign', 'sunjoong.TMalignWrap']

scripts = ['scripts/pyTMalign']


ext_name = 'sunjoong.TMalign'
sources = ['src/TMalign.pyf', 'src/TMalign_compares.f', 'src/TMalign.f']

include_dirs = []
library_dirs = []
libraries = []
define_macros = []
undef_macros = []
extra_objects = []
f2py_options = []


setup(name = name,
      version = version,
      author = author,
      author_email = author_email,
      description = description,
      long_description = long_description,
      url = url,
      license = license,
      classifiers = classifiers,
      platforms = 'All POSIX (Linux/BSD/UNIX-like OSes)',
      py_modules = py_modules,
      scripts = scripts,
      ext_modules = [Extension(name = ext_name,
                               sources = sources,
                               include_dirs = include_dirs,
                               library_dirs = library_dirs,
                               libraries = libraries,
                               define_macros = define_macros,
                               undef_macros = undef_macros,
                               extra_objects = extra_objects,
                               f2py_options = f2py_options)])

