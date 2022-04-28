# Copyright (c) 2006, Sunjoong LEE <sunjoong@gmail.com>
# Licensed under the Academic Free License version 3.0
# http://sunjoong.sourceforge.net/
__author__ = 'Sunjoong LEE <sunjoong@gmail.com>'
__date__ = '2006-07-08'
__version__ = '1.0.0'
__all__ = ['pyTMalign', 'pyTMalignComb', 'pyTMalignMatrix']

from FileUtils import openByFileExtension
from TMalign import tmalign, \
     backbone, pdb, length, options, filenames, d0, result, \
     init, alignrst, dpc, tm
from itertools import count, ifilter, imap, izip, dropwhile, takewhile
from struct import pack, unpack


_SLC = {'BCK':'X', 'GLY':'G', 'ALA':'A', 'SER':'S', 'CYS':'C',
        'VAL':'V', 'THR':'T', 'ILE':'I', 'PRO':'P', 'MET':'M',
        'ASP':'D', 'ASN':'N', 'LEU':'L', 'LYS':'K', 'GLU':'E',
        'GLN':'Q', 'ARG':'R', 'HIS':'H', 'PHE':'F', 'TYR':'Y',
        'TRP':'W', 'CYX':'C'}




def pyTMalign(structure, target, ave = False, sort = None,
              outname = None, L_fix = None, d0_min = None, d0_fix = None,
              filename1 = None, filename2 = None):
    """Example1:

    >>> structure = '153l.pdb'
    >>> target = '1dps.pdb'
    >>>
    >>> structure_size, target_size, normalized_size, aligned_length, \\
    ... rmsd, tm_score, seq_id, \\
    ... aligned_seq1, aligned_seq2, aligned_seq3 = \\
    ... pyTMalign(structure, target)
    >>>

    The return values - structure_size, target_size, normalized_size,
    aligned_length, rmsd, tm_score, seq_id, aligned_seq1, aligned_seq2,
    and aligned_seq3 - are variable result of 'TMalign' command.

    The arguments - structure and target - are pdb filenames
    or file pointers, etc.

    Example2:

    >>> size1, size2, normalized_size, aligned_length, \\
    ... rmsd, tm_score, seq_id, \\
    ... aligned_seq1, aligned_seq2, aligned_seq3 = \\
    ... pyTMalign(structure, target, ave = True)
    >>>

    Above example is equivalent to a '-a' option in 'TMalign' command.

    For a '-o TM.sup' option, it should be called by
    ... pyTMalign(structure, target, outname = 'TM.sup')

    If you want to use '-a' option and '-o TM.sup' option at once,
    ... pyTMalign(structure, target, ave = True, outname = 'TM.sup')

    Example3:

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

    Example4:

    >>> aligned = pyTMalign(structure, target, ave = True, sort = 'i')

    Above example use average length
    and return value will be sorted by sequence identity.

    Sort option is one of 'l', 'r', 's', 'i' and it's combinations.

    See also:
    To understand options details, read __TMalignWrap_getOpt() function code.
    """
    __pyTMalign_setOptions(ave, L_fix, d0_min, d0_fix)

    nseq2 = __pyTMalign_setTarget(target, filename2)
    only_one, structure, outname, filename1 = \
              __pyTMalign_fixList(structure, target, outname, filename1)
    aligned = []
    for (istructure, ioutname, ifilename1) in izip(
        structure, outname, filename1):
        __pyTMalign_setOutname(ioutname)
        nseq1 = __pyTMalign_setStructure(istructure, ifilename1)

        __pyTMalign_call_tmalign()

        nseq3, aligned_length, rmsd, tm_score, seq_id, \
               aligned_seq1, aligned_seq2, aligned_seq3 = \
               __pyTMalign_collectResults()

        aligned.append((nseq1, nseq2, nseq3, aligned_length,
                        rmsd, tm_score, seq_id,
                        aligned_seq1, aligned_seq2, aligned_seq3))

    if only_one:
        return aligned[0]
    else:
        if sort:
            __pyTMalign_sort(aligned, sort, offset = 3)
        return aligned




def pyTMalignComb(structure, ave = False, sort = None,
                  outname = None, L_fix = None, d0_min = None, d0_fix = None,
                  filename = None):
    """Example1:

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

    Example2:

    >>> aligned = pyTMalignComb(structure, ave = True, sort = 'i')

    Above example use average length
    and return value will be sorted by sequence identity.

    Sort option is one of 'l', 'r', 's', 'i' and it's combinations.

    See also:
    pyTMalign()
    """
    __pyTMalign_setOptions(ave, L_fix, d0_min, d0_fix)

    aligned = []
    for target, filename2, structure in \
        __pyTMalignComb_fixList(structure, outname, filename).__iter__():
        nseq2 = __pyTMalign_setTarget(target, filename2)

        for istructure, outname, filename1 in structure.__iter__():
            __pyTMalign_setOutname(outname)
            nseq1 = __pyTMalign_setStructure(istructure, filename1)
            __pyTMalign_call_tmalign()

            nseq3, aligned_length, rmsd, tm_score, seq_id, \
                   aligned_seq1, aligned_seq2, aligned_seq3 = \
                   __pyTMalign_collectResults()

            aligned.append((istructure, nseq1, target, nseq2,
                            nseq3, aligned_length, rmsd, tm_score, seq_id,
                            aligned_seq1, aligned_seq2, aligned_seq3))

    if sort:
        __pyTMalign_sort(aligned, sort, offset = 5)
    return aligned




def pyTMalignMatrix(structure, ave = False,
                    outname = None, L_fix = None, d0_min = None, d0_fix = None,
                    filename = None):
    """Example:

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

    Caution:
    pyTMalignMatrix() does not take sort argument.

    See also:
    pyTMalign()
    """
    __pyTMalign_setOptions(ave, L_fix, d0_min, d0_fix)

    itarget = ''
    for istructure, target, ioutname, filename1, filename2 in \
            __pyTMalignMatrix_fixList(structure, outname, filename):
        if target != itarget:
            nseq2 = __pyTMalign_setTarget(target, filename2)
            itarget = target

        __pyTMalign_setOutname(ioutname)
        nseq1 = __pyTMalign_setStructure(istructure, filename1)

        __pyTMalign_call_tmalign()

        nseq3, aligned_length, rmsd, tm_score, seq_id, \
               aligned_seq1, aligned_seq2, aligned_seq3 = \
               __pyTMalign_collectResults()

        yield (istructure, nseq1, target, nseq2, nseq3,
               aligned_length, rmsd, tm_score, seq_id,
               aligned_seq1, aligned_seq2, aligned_seq3)
    raise StopIteration




def __pyTMalign_setOptions(ave, L_fix, d0_min, d0_fix):
    if ave:
        options.m_ave.data.__setslice__(0, 4, pack('i', 1))
    else:
        options.m_ave.data.__setslice__(0, 4, pack('i', -1))
    if L_fix:
        options.m_fix.data.__setslice__(0, 4, pack('i', 1))
        options.l_fix.data.__setslice__(0, 4, pack('i', L_fix))
    else:
        options.m_fix.data.__setslice__(0, 4, pack('i', -1))
        options.l_fix.data.__setslice__(0, 4, pack('i', 0))
    if d0_min:
        options.m_d0_min.data.__setslice__(0, 4, pack('i', 1))
        options.d0_min_input.data.__setslice__(0, 4, pack('f', d0_min))
    else:
        options.m_d0_min.data.__setslice__(0, 4, pack('i', -1))
        options.d0_min_input.data.__setslice__(0, 4, pack('f', 0.0))
    if d0_fix:
        options.m_d0.data.__setslice__(0, 4, pack('i', 1))
        options.d0_fix.data.__setslice__(0, 4, pack('f', d0_fix))
    else:
        options.m_d0.data.__setslice__(0, 4, pack('i', -1))
        options.d0_fix.data.__setslice__(0, 4, pack('f', 0.0))


def __pyTMalign_setOutname(outname):
    if outname:
        options.m_out.data.__setslice__(0, 4, pack('i', 1))
        filenames.outname.data.__setslice__(0, len(outname), outname)
    else:
        options.m_out.data.__setslice__(0, 4, pack('i', -1))




def __pyTMalign_setTarget(target, filename):
    filename, nseq =  __pyTMalign_readPdb(target, filename,
                                          backbone.xa, 1,
                                          pdb.ss2, pdb.mm2, pdb.seq2)
    filenames.pdb2.data.__setslice__(0, len(filename), filename)
    length.nseq2.data.__setslice__(0, 4, pack('i', nseq))
    return nseq


def __pyTMalign_setStructure(structure, filename):
    filename, nseq =  __pyTMalign_readPdb(structure, filename,
                                          backbone.xa, 0,
                                          pdb.ss1, pdb.mm1, pdb.seq1)
    filenames.pdb1.data.__setslice__(0, len(filename), filename)
    length.nseq1.data.__setslice__(0, 4, pack('i', nseq))
    return nseq


def __pyTMalign_readPdb(pdb, filename, xa, xa_index, ss, mm, seq):
    if type(pdb) is type('') and pdb:
        fp = openByFileExtension(pdb)
        nseq = __pyTMalign_readPdb_fromFile(fp, xa, xa_index, ss, mm, seq)
        fp.close()
        return (pdb, nseq)

    elif hasattr(pdb, 'readlines'):
        nseq = __pyTMalign_readPdb_fromFile(pdb, xa, xa_index, ss, mm, seq)
        if hasattr(pdb, 'name') and pdb.name:
            return (pdb.name, nseq)
        elif filename:
            return (filename, nseq)
        else:
            return ('', nseq)

    elif hasattr(pdb, '_is_a_typed_list'):
        name, nseq = __pyTMalign_readPdb_fromMol(pdb, xa, xa_index,
                                                 ss, mm, seq)
        if filename:
            return (filename, nseq)
        else:
            return (name, nseq)

    else:
        raise ValueError, 'Unknown type: %s.' % pdb


def __pyTMalign_readPdb_fromFile(fp, xa, xa_index, ss, mm, seq):
    fp.tell() or fp.seek(0)
    return len(map(lambda (x, y): (

        ss.__setitem__(y, x[0]),
        seq.data.__setitem__(y + 1, _SLC[x[0]]), # caution! seq begins from 0.
        mm.__setitem__(y, x[1]),
        xa[0, y].__setitem__(xa_index, x[2]),
        xa[1, y].__setitem__(xa_index, x[3]),
        xa[2, y].__setitem__(xa_index, x[4])),

                   izip(
        imap(lambda x: (x[17:20], eval(x[23:26]),
                        eval(x[31:38]), eval(x[39:46]), eval(x[47:54])),

             ifilter(lambda x: (x[12:15].strip() == 'CA'
                                and (x[16] == ' ' or x[16] == 'A')
                                and x[:5] == 'ATOM '),

                     takewhile(lambda x: x[:3] != 'TER',
                               dropwhile(lambda x: x[:5] != 'ATOM ',
                                         fp.readlines().__iter__())))),

        count())))


def __pyTMalign_readPdb_fromMol(pdb, xa, xa_index, ss, mm, seq):
    if hasattr(pdb, '_is_a_residue'):
        raise ValueError, '%s is a residue' % pdb
    elif hasattr(pdb, '_is_a_chain'):
        name = (pdb.name == ' ' and 'chain_'
                or 'chain%s' % pdb.name.upper())
    else:
        name = pdb.name[:100]
        pdb = pdb[0]

    size = len(map(lambda (x, y): (

        ss.__setitem__(y, pdb[y].name),
        seq.data.__setitem__(y + 1, _SLC[pdb[y].name]),
        mm.__setitem__(y, x[0]),
        xa[0, y].__setitem__(xa_index, x[1]),
        xa[1, y].__setitem__(xa_index, x[2]),
        xa[2, y].__setitem__(xa_index, x[3])),

                   izip(
        imap(lambda x: (x.idx, x.x, x.y, x.z),
             pdb.atoms_with_name('CA')),
        count())))

    return (name, size)




def __pyTMalign_fixList(structure, target, outname, filename):
    if type(structure) == type(''):
        if outname and type(outname) == type([]):
            outname = outname[0]
        elif outname is not None:
            raise ValueError, 'outname %s' % outname
        if filename and type(filename) == type([]):
            filename = filename[0]
        elif filename is not None:
            raise ValueError, 'filename %s' % filename
        return (True, [structure], [outname], [filename])
    elif type(structure) == type([]):
        if outname and type(outname) == type(''):
            base_name = outname
            outname = map(lambda x: '%s_%s_%s' % (base_name, x, target),
                          structure.__iter__())
        elif outname is None:
            outname = [None] * len(structure)
        elif type(outname) != type([]) or len(outname) != len(structure):
            raise ValueError, 'outname %s' % outname
        if filename and type(filename) == type(''):
            base_name = filename
            filename = map(lambda x: '%s_%s' % (base_name, x),
                           xrange(len(structure)))
        elif filename is None:
            filename = [None] * len(structure)
        elif type(filename) != type([]) or len(filename) != len(structure):
            raise ValueError, 'filename %s' % filename
        return (False, structure, outname, filename)
    else:
        raise ValueError, 'structure %s' % structure


def __pyTMalignComb_fixList(structure, outname, filename):
    if len(structure) < 2:
        raise ValueError, 'structure size error: %s' % structure

    if outname and type(outname) == type(''):
        getOutname = (lambda x, y: '%s_%s_%s' % (outname, x, y))
    elif outname is None:
        getOutname = (lambda x, y: None)
    else:
        raise ValueError, 'outname should be string or None'

    if filename and type(filename) == type(''):
        getFilename = (lambda x: '%s_%s' % (filename, x))
    elif filename is None:
        getFilename = (lambda x: None)
    else:
        raise ValueError, 'filename should be string or None'

    return [(structure[x:][0], getFilename(structure[x:][0]),
             [(y, getOutname(y, structure[x:][0]), getFilename(y))
              for y in structure[:x]])
            for x in xrange(len(structure) - 1, 0, -1)]


def __pyTMalignMatrix_fixList(structure, outname, filename):
    if len(structure) < 2:
        raise ValueError, 'structure size error: %s' % structure

    if outname and type(outname) == type(''):
        getOutname = (lambda x, y: '%s_%s_%s' % (outname, x, y))
    elif outname is None:
        getOutname = (lambda x, y: None)
    else:
        raise ValueError, 'outname should be string or None'

    if filename and type(filename) == type(''):
        getFilename = (lambda x: '%s_%s' % (filename, x))
    elif filename is None:
        getFilename = (lambda x: None)
    else:
        raise ValueError, 'filename should be string or None'

    for target in structure.__iter__():
        for istructure in structure.__iter__():
            yield (istructure, target, getOutname(istructure, target),
                   getFilename(istructure), getFilename(target))
    raise StopIteration




def __pyTMalign_call_tmalign():
    init.invmap_i.__setslice__(0, 3000, 0)
    alignrst.invmap0.__setslice__(0, 3000, 0)
    dpc.invmap.__setslice__(0, 3000, 0)
    map(lambda x: x.__setslice__(0, 3000, 0), dpc.score)
    dpc.gap_open.data.__setslice__(0, 4, pack('f', 0.0))
    tm.tm.data.__setslice__(0, 4, pack('f', 0))
    tm.tmmax.data.__setslice__(0, 4, pack('f', 0))

    result.aseq1.data.__setslice__(0, 6000, pack('c', '\0') * 6000)
    result.aseq2.data.__setslice__(0, 6000, pack('c', '\0') * 6000)
    result.aseq3.data.__setslice__(0, 6000, pack('c', '\0') * 6000)

    tmalign()


def __pyTMalign_collectResults():
    nseq3 = int(unpack('f', d0.anseq.data)[0])
    aligned_length = unpack('i', result.n8_al.data)[0]
    rmsd = unpack('f', result.rmsd8_al.data)[0]
    tm_score = unpack('f', result.tm8.data)[0]
    seq_id = unpack('f', result.seq_id.data)[0]
    aligned_seq1 = result.aseq1.tostring().strip('\x00')
    aligned_seq2 = result.aseq2.tostring().strip('\x00')
    aligned_seq3 = result.aseq3.tostring().strip('\x00')

    return (nseq3, aligned_length, rmsd, tm_score, seq_id,
            aligned_seq1, aligned_seq2, aligned_seq3)




def __pyTMalign_sort(aligned, sort, offset = 3):
    sort = list(sort)
    sort.reverse()
    for i in sort:
        if i == 'l':
            aligned.sort(lambda x, y: cmp(x[offset], y[offset]))
        elif i == 'r':
            aligned.sort(lambda x, y: cmp(x[offset + 1], y[offset + 1]))
        elif i == 's':
            aligned.sort(lambda x, y: cmp(x[offset + 2], y[offset + 2]))
        elif i == 'i':
            aligned.sort(lambda x, y: cmp(x[offset + 3], y[offset + 3]))
        else:
            raise ValueError, 'sort %s for %s' % (sort, i)

