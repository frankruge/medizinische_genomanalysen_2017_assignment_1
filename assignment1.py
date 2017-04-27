#! /usr/bin/env python3
#from __future__ import division
import mysql.connector
import os
import re
import pysam
import pybedtools as pbt
import numpy as np
import subprocess
import string
import decimal

os.chdir("/home/frank/Desktop/FH_Semester2/medizinischeGenomanalysepapinger/medizinische_genomanalysen_2017_assignment_1/")
__author__ = 'Frank Ruge'
pick_name='INS'
#pick_name=raw_input("pick gene and press enter:")

annotation_file=''

class Assignment1:
    
    def __init__(self, name):
        ## Your gene of interest
        self.name = name
        self.genes = open("TEST.TXT", 'r')
        print("fetching data for: " + self.name)

        self.bamfile = os.path.join(os.getcwd(), "HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.sorted.bam")
        #pysam.sort(self.bamfile, "HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.sorted.bam")
        #self.sorted_bamfile = os.path.join(os.getcwd(), "HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.sorted.bam")
        #pysam.index(self.sorted_bamfile)
        self.bf = pysam.AlignmentFile(self.bamfile, 'rb')



    
    def fetch_gene_coordinates(self, genome_reference, file_name):
        #if the GENE was searched before, a file GENE.txt was generated in os.getcwd(). it contains the coordinates. (TEST.TXT if you cloned the repository)
        # the entire UCSC part of this script is omitted if the file exists, thereby saving time. in order to run the
        # UCSC part of this script, GENE.txt or TEST.TXT has to be removed.
        if str(self.name + ".txt") in os.listdir(os.getcwd()):
            mygene=[]
            coords= open(str(self.name +  ".txt"), 'r')
            a=coords.readlines()
            for i in range(len(a)):
                a[i]=a[i].strip('\n')
            coords.close()
            mygene.append(a)
            if mygene[0]:
                #print(self.name+" has these coordinates: "+str(mygene[0][2].strip('\n'))+'\t'+str(mygene[0][3].strip('\n'))+'\t'+str(mygene[0][4].strip('\n'))+"\nYou have previously collected the coordinates for this gene\n")
                return mygene
            #else:
             #   pass

        print("Connecting to UCSC to fetch data")
        #print(self.name)

        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)

        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        ## Execute query
        cursor.execute(query)
        ## Write to file
        mygene=[]
        if mygene:
            print("variable is not null")
        fh=open(file_name, 'w')
        for row in cursor:
            if self.name in row:
                if mygene:
                    if mygene[-1][3] == row[3] and mygene[-1][4] == row[4] and mygene[-1][2] == row[2]:
                        print(self.name+" multiple entries. omitting entries with same coordinates")
                        pass
                else:
                    mygene.append(row)

                fh.write(str(row) + "\n")
            fh.write(str(row) + "\n")
        fh.close()

        ## Close cursor & connection
        cursor.close()
        cnx.close()
        if len(mygene)>1:
            a=input('which list do you choose (1,2,etc.)?\n'+str(mygene)) #user chooses from multiple entries
            a=a-1
        else:
            a=0
        try:    # test if mygene[] has at least one entry. if not the gene name given by user doesnt exist in the ucsc list
            coordinates=[mygene[a][2],mygene[a][3],mygene[a][4]]
        except:
            return "your gene seems not to be listed. doublecheck name, capitalization, spacing, species, etc!"
        print("Done fetching data")
        print(mygene[a][0]+ " coordinates are: "+str(mygene[a][2])+'\t'+str(mygene[a][3])+'\t'+str(mygene[a][4]))
        print("Created file containing all entries retrieved from ucsc: " + file_name)
        cfile= open(str(self.name + ".txt"), "w")
        for entry in mygene[a]:
            #cfile.write(mygene[a])
            cfile.write(str(entry)+'\n')
        #cfile.write(str(str(coordinates[2])+'\t'+str(coordinates[3])+'\t'+str(coordinates[4])))
        cfile.close()
        return mygene

    def get_sam_header(self):
        self.bf = pysam.AlignmentFile(self.bamfile, 'rb')
        #header = pysam.view(self.bamfile, '-H') #is like samtools view -H file.bam
        header = self.bf.header #returns a two level dictionary with information
        #for list in header.items():
        #    print(list[0]+'\t'+str([i for i in list[1]]))
        self.bf.close()
        for i in header:
            print("header: "+str(header[i]))
        return(header)

    def get_properly_paired_reads_of_gene(self, outname, coordinates):
        self.bf = pysam.AlignmentFile(self.bamfile, 'rb')
        pairedreads= pysam.AlignmentFile(outname, 'wb', template=self.bf)
        count=0

        for read in self.bf.fetch(reference=str(coordinates[0][2]).strip('chr'), start=int(coordinates[0][3]), end=int(coordinates[0][4])):
            if read.is_proper_pair:
                pairedreads.write(read)
                count+=1
        pairedreads.close()
        #print("\nNumber of proper pairs at " + self.name + " coordinates" + str(coordinates)+":\n" + str(count) + "\n\n")
        print("get_properly_paired_reads_of_gene: " + str(count))
        self.bf.close()
        return

        
    def get_gene_reads_with_indels(self, coordinates):
        self.bf = pysam.AlignmentFile(self.bamfile, 'rb')
        indels=[]
        count=0
        for read in self.bf.fetch(reference=str(coordinates[0][2]).strip('chr'), start=int(coordinates[0][3]), end=int(coordinates[0][4])):
            if not read.is_unmapped:
                cigar_line = read.cigartuples

                for entry in cigar_line:
                    if len(cigar_line) > 1:
                        #indels.append(read)
                        count+=1
        self.bf.close()
        #print("Number of indels found at the region " + self.name + " coordinates" + str(coordinates) +": \n" + str(len(indels)))
        print("get_gene_reads_with_indels: " + str(count))
        return #indels


    def calculate_total_average_coverage(self, readlength):
        #number of reads * length of each read / size of the genome
        self.bf = pysam.AlignmentFile(self.bamfile, 'rb')
        counter=0
        genome=0
        #cov=self.bf.count_coverage() #pysam function has a bug. it is documented but has never been resolved
        for read in self.bf:
            if not read.is_unmapped:
                counter+=1
        #genome size = sum length of chromosomes
        for entry in self.bf.lengths:
            genome = genome + entry
        coverage = ((counter*readlength)/genome)
        print("calculate_total_average_coverage "+str(round(coverage, 3)))
        self.bf.close()
        return

    def calculate_gene_average_coverage(self, coordinates, readlength):
        self.bf = pysam.AlignmentFile(self.bamfile, 'rb')
        subset=self.bf.fetch(reference=str(coordinates[0][2]).strip('chr'), start=int(coordinates[0][3]), end=int(coordinates[0][4]))
        counter = 0
        gene =  int(coordinates[0][4]) - int(coordinates[0][3])#length og gene
        # cov=self.bf.count_coverage() #pysam function has a bug. it is documented but has never been resolved
        for read in subset:
            if not read.is_unmapped:
                counter += 1
        coverage = ((counter * readlength) / gene)
        print("calculate_gene_average_coverage " + str(round(coverage, 3)))
        self.bf.close()
        return

    def get_number_mapped_reads(self):
        self.bf = pysam.AlignmentFile(self.bamfile, 'rb')
        count=0
        for read in self.bf:
            if not read.is_unmapped:
                count+=1
        print("get_number_mapped_reads: "+str(count))
        self.bf.close()
        return (count)
        
    def get_gene(self, symbol):
        gene_line=[]
        sym="'"+symbol+"'"
        self.genes = open("TEST.TXT", 'r')
        for line in self.genes:
            if sym in line:
                a = re.split(", ", line.strip('()').replace("'", ""))
                if gene_line:                         #if gene_line has an entry
                    if a[3] == gene_line[-1][3] and a[4] == gene_line[-1][4]:        #check if new entry has the same coordinates as previous entry
                        pass
                else:
                    gene_line.append(a)
        #print("get_gene: "+str(gene_line))
        self.genes.close()
        return(gene_line)

    def get_gene_symbol(self, symbol):
        #self.genes = open("TEST.TXT", 'r')
        gl = assignment1.get_gene(symbol)
        print("get gene symbol: "+gl[0][0])
        #print("get gene symbol: " + gl[0][5])
        #print("get gene symbol: "+gl[0][6])
        #print("get gene symbol: "+gl[0][7])
        #self.genes.close()
        return(gl[0][0])

        
    def get_region_of_gene(self, symbol):
        #self.genes = open("TEST.TXT", 'r')
        gl = assignment1.get_gene(symbol)
        coords=[gl[0][2],gl[0][3], gl[0][4]]
        #coords = gl
        print("get_region_of_gene" + str(coords))
        #self.genes.close()
        return coords
        
    def get_number_of_exons(self, symbol):

        gl = assignment1.get_gene(symbol)
        exons = gl[0][6]
        print("number of exons: " + str(exons))
        return exons
    
    def print_summary(self):
        coordinates=tuple(assignment1.fetch_gene_coordinates("hg19", "human_genes.table")); #print(coordinates) ##
        header = assignment1.get_sam_header();#print(header)
        #getting paired reads takes 2 min for example file
        proper_paired=assignment1.get_properly_paired_reads_of_gene(outname="paired.bam", coordinates=coordinates)
        indels=assignment1.get_gene_reads_with_indels(coordinates=coordinates)
        cov=assignment1.calculate_total_average_coverage(readlength=100);#print(cov)
        assignment1.calculate_gene_average_coverage(coordinates=coordinates, readlength=100)
        mapped=assignment1.get_number_mapped_reads();#print(mapped)
        gene_symbol=assignment1.get_gene_symbol(symbol=self.name)#"INS");#print(gene_symbol)
        region=assignment1.get_region_of_gene(symbol=self.name);#print(region)
        exons = assignment1.get_number_of_exons(symbol=self.name);#print(exons)
    
        
if __name__ == '__main__':
    print("Assignment 1")
    print(__author__)

    assignment1 = Assignment1(name=pick_name)

    assignment1.print_summary()