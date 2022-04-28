# -*- coding: utf-8 -*-
from test import test_support
from unittest import TestCase, TestSuite, makeSuite

from sunjoong.pyTMalign import pyTMalign, pyTMalignComb, pyTMalignMatrix

from cPickle import load as pickle_load
from itertools import ifilter


class pyTMalign_Test(TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_1(self):
        """수동으로 검사했던 결과와 비교"""
        pdb = ['153l', '1dps', '1e70', '1qp8', '1xdw', '2fcf']
        for structure, target in ifilter(lambda (x, y): x != y,
                                         [(x, y)
                                          for x in pdb
                                          for y in pdb].__iter__()):
            fp = open('data/%s_%s.dump' % (structure, target), 'r')
            expected = pickle_load(fp)
            fp.close()

            result = pyTMalign('data/%s.pdb' % structure,
                               'data/%s.pdb' % target)

            self.assertEquals(len(expected), len(result))
            map(lambda x, y: self.assertEquals(x, y),
                expected.__iter__(), result.__iter__())

            fp = open('data/%s_%s_ave.dump' % (structure, target), 'r')
            expected = pickle_load(fp)
            fp.close()

            result = pyTMalign('data/%s.pdb' % structure,
                               'data/%s.pdb' % target, ave = True)

            self.assertEquals(len(expected), len(result))
            map(lambda x, y: self.assertEquals(x, y),
                expected.__iter__(), result.__iter__())


    def test_2(self):
        """structure가 리스트 일때"""
        structure = ['data/153l.pdb', 'data/1dps.pdb', 'data/1e70.pdb',
                     'data/1qp8.pdb', 'data/1xdw.pdb']
        target = 'data/2fcf.pdb'

        expected = []
        for i in ['data/153l_2fcf.dump', 'data/1dps_2fcf.dump',
                  'data/1e70_2fcf.dump', 'data/1qp8_2fcf.dump',
                  'data/1xdw_2fcf.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target)

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())

        expected = []
        for i in ['data/153l_2fcf_ave.dump', 'data/1dps_2fcf_ave.dump',
                  'data/1e70_2fcf_ave.dump', 'data/1qp8_2fcf_ave.dump',
                  'data/1xdw_2fcf_ave.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target, ave = True)

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())


    def test_sort_1(self):
        """소트 시험1"""
        structure = ['data/153l.pdb', 'data/1dps.pdb', 'data/1e70.pdb',
                     'data/1qp8.pdb', 'data/1xdw.pdb']
        target = 'data/2fcf.pdb'

        expected = []
        for i in ['data/1qp8_2fcf.dump', 'data/1dps_2fcf.dump',
                  'data/153l_2fcf.dump', 'data/1e70_2fcf.dump',
                  'data/1xdw_2fcf.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target, sort = 'i')

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())

        expected = []
        for i in ['data/1qp8_2fcf_ave.dump', 'data/1dps_2fcf_ave.dump',
                  'data/153l_2fcf_ave.dump', 'data/1e70_2fcf_ave.dump',
                  'data/1xdw_2fcf_ave.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target, ave = True, sort = 'i')

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())


    def test_sort_2(self):
        """소트 시험2"""
        structure = ['data/153l.pdb', 'data/1dps.pdb', 'data/1e70.pdb',
                     'data/1qp8.pdb', 'data/1xdw.pdb']
        target = 'data/2fcf.pdb'

        expected = []
        for i in ['data/1dps_2fcf.dump', 'data/153l_2fcf.dump',
                  'data/1e70_2fcf.dump', 'data/1xdw_2fcf.dump',
                  'data/1qp8_2fcf.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target, sort = 's')

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())

        expected = []
        for i in ['data/1e70_2fcf_ave.dump', 'data/1dps_2fcf_ave.dump',
                  'data/1xdw_2fcf_ave.dump', 'data/1qp8_2fcf_ave.dump',
                  'data/153l_2fcf_ave.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target, ave = True, sort = 's')

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())


    def test_sort_3(self):
        """소트 시험3"""
        structure = ['data/153l.pdb', 'data/1dps.pdb', 'data/1e70.pdb',
                     'data/1qp8.pdb', 'data/1xdw.pdb']
        target = 'data/2fcf.pdb'

        expected = []
        for i in ['data/1qp8_2fcf.dump', 'data/153l_2fcf.dump',
                  'data/1dps_2fcf.dump', 'data/1xdw_2fcf.dump',
                  'data/1e70_2fcf.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target, sort = 'r')

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())

        expected = []
        for i in ['data/1qp8_2fcf_ave.dump', 'data/153l_2fcf_ave.dump',
                  'data/1dps_2fcf_ave.dump', 'data/1xdw_2fcf_ave.dump',
                  'data/1e70_2fcf_ave.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target, ave = True, sort = 'r')

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())


    def test_sort_4(self):
        """소트 시험4"""
        structure = ['data/153l.pdb', 'data/1dps.pdb', 'data/1e70.pdb',
                     'data/1qp8.pdb', 'data/1xdw.pdb']
        target = 'data/2fcf.pdb'

        expected = []
        for i in ['data/1dps_2fcf.dump', 'data/153l_2fcf.dump',
                  'data/1qp8_2fcf.dump', 'data/1e70_2fcf.dump',
                  'data/1xdw_2fcf.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target, sort = 'l')

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())

        expected = []
        for i in ['data/1dps_2fcf_ave.dump', 'data/153l_2fcf_ave.dump',
                  'data/1qp8_2fcf_ave.dump', 'data/1e70_2fcf_ave.dump',
                  'data/1xdw_2fcf_ave.dump'].__iter__():
            fp = open(i, 'r')
            expected.append(pickle_load(fp))
            fp.close()

        result = pyTMalign(structure, target, ave = True, sort = 'l')

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())




class pyTMalignComb_Test(TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_1(self):
        """수동으로 검사했던 결과와 비교"""
        structure = ['data/153l.pdb', 'data/1dps.pdb', 'data/1e70.pdb',
                     'data/1qp8.pdb', 'data/1xdw.pdb', 'data/2fcf.pdb']

        expected = []
        for i, j in [('153l', '2fcf'), ('1dps', '2fcf'), ('1e70', '2fcf'),
                     ('1qp8', '2fcf'), ('1xdw', '2fcf'), ('153l', '1xdw'),
                     ('1dps', '1xdw'), ('1e70', '1xdw'), ('1qp8', '1xdw'),
                     ('153l', '1qp8'), ('1dps', '1qp8'), ('1e70', '1qp8'),
                     ('153l', '1e70'), ('1dps', '1e70'), ('153l', '1dps')
                     ].__iter__():
            istructure = 'data/%s.pdb' % i
            target = 'data/%s.pdb' % j
            fp = open('data/%s_%s.dump' % (i, j), 'r')
            aligned = list(pickle_load(fp))
            fp.close()
            aligned.insert(0, istructure)
            aligned.insert(2, target)
            expected.append(aligned)

        result = pyTMalignComb(structure)

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())

        expected = []
        for i, j in [('153l', '2fcf'), ('1dps', '2fcf'), ('1e70', '2fcf'),
                     ('1qp8', '2fcf'), ('1xdw', '2fcf'), ('153l', '1xdw'),
                     ('1dps', '1xdw'), ('1e70', '1xdw'), ('1qp8', '1xdw'),
                     ('153l', '1qp8'), ('1dps', '1qp8'), ('1e70', '1qp8'),
                     ('153l', '1e70'), ('1dps', '1e70'), ('153l', '1dps')
                     ].__iter__():
            istructure = 'data/%s.pdb' % i
            target = 'data/%s.pdb' % j
            fp = open('data/%s_%s_ave.dump' % (i, j), 'r')
            aligned = list(pickle_load(fp))
            fp.close()
            aligned.insert(0, istructure)
            aligned.insert(2, target)
            expected.append(aligned)

        result = pyTMalignComb(structure, ave = True)

        self.assertEquals(len(expected), len(result))
        map(lambda x, y: map(lambda a, b: self.assertEquals(a, b),
                             x.__iter__(), y.__iter__()),
            expected.__iter__(), result.__iter__())




class pyTMalignMatrix_Test(TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_1(self):
        """수동으로 검사했던 결과와 비교"""
        structure = ['data/153l.pdb', 'data/1dps.pdb', 'data/1e70.pdb',
                     'data/1qp8.pdb', 'data/1xdw.pdb', 'data/2fcf.pdb']

        kount = 0
        for istructure, structure_size, target, target_size, \
                normalized_size, aligned_length, rmsd, tm_score, seq_id, \
                aligned_seq1, aligned_seq2, aligned_seq3 \
                in pyTMalignMatrix(structure):
            kount += 1
            if istructure != target:
                fp = open('data/%s_%s.dump' % (istructure[5:-4],
                                               target[5:-4]), 'r')
                expected = pickle_load(fp)
                fp.close()

                result = (structure_size, target_size, normalized_size,
                          aligned_length, rmsd, tm_score, seq_id,
                          aligned_seq1, aligned_seq2, aligned_seq3)

                map(lambda x, y: self.assertEquals(x, y),
                    expected.__iter__(), result.__iter__())
        self.assertEquals(len(structure) ** 2, kount)

        kount = 0
        for istructure, structure_size, target, target_size, \
                normalized_size, aligned_length, rmsd, tm_score, seq_id, \
                aligned_seq1, aligned_seq2, aligned_seq3 \
                in pyTMalignMatrix(structure, ave = True):
            kount += 1
            if istructure != target:
                fp = open('data/%s_%s_ave.dump' % (istructure[5:-4],
                                                   target[5:-4]), 'r')
                expected = pickle_load(fp)
                fp.close()

                result = (structure_size, target_size, normalized_size,
                          aligned_length, rmsd, tm_score, seq_id,
                          aligned_seq1, aligned_seq2, aligned_seq3)

                map(lambda x, y: self.assertEquals(x, y),
                    expected.__iter__(), result.__iter__())
        self.assertEquals(len(structure) ** 2, kount)




def TestMain():
    suite = TestSuite()
    suite.addTest(makeSuite(pyTMalign_Test))
    suite.addTest(makeSuite(pyTMalignComb_Test))
    suite.addTest(makeSuite(pyTMalignMatrix_Test))
    test_support.run_suite(suite)




if __name__ == '__main__':
    TestMain()

