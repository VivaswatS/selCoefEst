#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

#*************************************************************************
#    > File Name: Ms2Vcf.py
#    > Author: xlzh
#    > Mail: xiaolongzhang2015@163.com 
#    > Created Time: 2019年08月25日 星期日 15时46分09秒
#*************************************************************************

import pdb
import sys
import numpy as np

class MS(object):
    def __init__(self, ms_file, sam_num, num_sites):
        self.ms_fp   = open(ms_file, 'r') # file handle of ms file
        self.pos     = self.__get_pos(num_sites)   # get the position list of ms file
        self.pos_num = len(self.pos)      # position number
        self.sam_num = sam_num            # haplotypes number
        self.geno    = np.zeros((self.sam_num, self.pos_num), dtype='int8')
        
    def __get_pos(self, num_sites):
        ''' get the position list of ms/msms file
        '''
        pos_list = []

        while (1):
            line = self.ms_fp.readline()
            if line.startswith('positions'):
                pos_list = line.split(':')[1].split()
                pos_list = [int(float(p)*num_sites) for p in pos_list]
                break

        return pos_list 


def main():
    args = sys.argv

    if len(args) != 4:
        sys.stderr.write('Usage: python %s <in.ms> <haps_num>\n' %args[0])
        sys.exit(-1)

    ms = MS(args[1], int(args[2]), int(args[3]))

    idx = 0
    # read the ms/msms file to get the haplotype
    for line in ms.ms_fp:
        if line.rstrip() == '': continue
        g_list = list(line.rstrip())
        for i in range(ms.pos_num): ms.geno[idx][i] = g_list[i]
        idx += 1

    # output vcf file
    out_fp = open('infiles/' + args[1].split('.ms')[0].split('/')[-1] + '.vcf', 'w')

    out_fp.write('##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">\n##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##contig=<ID=1,length=2147483647>\n')

    header = '\t'.join([str(i+1) for i in range(ms.sam_num//2)])
    out_fp.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' %header)

    # construct the genotype with nearby two haplotypes
    for p_idx in range(ms.pos_num):
        sam_geno = []
        for s_idx in range(0, ms.sam_num, 2):
            sam_geno.append(str(ms.geno[s_idx][p_idx])+'|'+str(ms.geno[s_idx+1][p_idx]))

        out_fp.write('1\t%d\tSNP%d\tA\tT\t.\tPASS\t.\tGT\t%s\n' %(ms.pos[p_idx],p_idx,'\t'.join(sam_geno)))

    out_fp.close()

if __name__ == '__main__':
    main()
