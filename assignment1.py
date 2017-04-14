import mysql.connector
import os
import os
import re
import pysam
import pybedtools as pbt
import numpy as np
import subprocess
import string

os.chdir("/home/frank/Desktop/FH_Semester2/medizinischeGenomanalysepapinger/medizinische_genomanalysen_2017_assignment_1/")
__author__ = 'Frank Ruge'
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
        fh=open(file_name, 'w')
        for row in cursor:
            #if self.name in row:
            #    fh.write(str(row) + "\n")
            fh.write(str(row) + "\n")


        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print("Done fetching data")
        return(fh)
    def get_sam_header(self):
#        header = pysam.view(self.bamfile, '-H') #is like samtools view -H file.bam
        header = self.bf.header #returns a two level dictionary with information
        for list in header.items():
            print(list[0]+'\t'+str([i for i in list[1]]))
        return(header)

    def get_properly_paired_reads_of_gene(self, outname):
        import pysam
        pairedreads= pysam.AlignmentFile(outname, 'wb', template=self.bf)
        for read in self.bf.fetch():
            if read.is_paired:
                pairedreads.write(read)
        pairedreads.close()

        
    def get_gene_reads_with_indels(self):
        indels=[]
        for read in self.bf:
            if not read.is_unmapped:
                cigar_line = read.cigartuples

                for entry in cigar_line:
                    if len(cigar_line) > 1:
                        indels.append(read)
                        #print(cigar_line)
        return(indels)
        print("todo")
        
    def calculate_total_average_coverage(self, readlength):
        #number of reads * length of each read / size of the genome
        count=0
        genome=0
        #cov=self.bf.count_coverage() #pysam function has a bug. it is documented but has never been resolved
        for read in self.bf:
            if not read.is_unmapped:
                count+=1

        #genome size = sum length of chromosomes
        for entry in self.bf.lengths:
            genome=genome+entry
        return((count*readlength)/genome)

    def calculate_gene_average_coverage(self):
        print("todo")
        
    def get_number_mapped_reads(self):
        count=0
        for read in self.bf:
            if not read.is_unmapped:
                count+=1
        return(count)
        
    def get_gene(self, symbol):
        gene_line=[]
        sym="'"+symbol+"'"
        genes=self.genes
        for line in genes:
            if sym in line:
                a = re.split(", ", line.strip('()').replace("'", ""))
                if gene_line:                         #if gene_line has an entry
                    if a[3] == gene_line[-1][3] and a[4] == gene_line[-1][4]:        #check if new entry has the same coordinates as previous entry
                        pass
                else:
                    gene_line.append(a)
        return(gene_line)

    def get_gene_symbol(self, symbol):
        gl=assignment1.get_gene(symbol)
        return(gl[0][1])

        
    def get_region_of_gene(self, symbol):
        gl = assignment1.get_gene(symbol)
        coords=[gl[0][3], gl[0][4]]
        return (coords)
        
    def get_number_of_exons(self, symbol):
        gl = assignment1.get_gene(symbol)
        exons = len(gl[0][7].strip(',').split(','))
        return exons
        print("ads")
    
    def print_summary(self):
        fh=assignment1.fetch_gene_coordinates("hg19", "human_genes.table") ##
        # header=assignment1.get_sam_header();print(header)
        # getting paired reads takes 2 min for example file
        # proper_paired=assignment1.get_properly_paired_reads_of_gene("paired.bam")

        # indel count is 2 787 368
        # indels=assignment1.get_gene_reads_with_indels(); print(len(indels))
        # cov=assignment1.calculate_total_average_coverage(readlength=100);print(cov)
        # mapped=assignment1.get_number_mapped_reads();print(mapped)

        # gene_symbol=assignment1.get_gene_symbol(symbol="INS");print(gene_symbol)
        # region=assignment1.get_region_of_gene("INS");print(region)
        exons = assignment1.get_number_of_exons("INS");print(exons)
    
        
if __name__ == '__main__':
    print("Assignment 1")
    print(__author__)
    assignment1 = Assignment1(name="INS")
    #assignment1.print_summary()
    
    #assignment1.fetch_gene_coordinates("hg19", "human_genes.table") ##
    #header=assignment1.get_sam_header();print(header)

    #getting paired reads takes 2 min for example file
    #proper_paired=assignment1.get_properly_paired_reads_of_gene("paired.bam")

    #indel count is 2 787 368
    #indels=assignment1.get_gene_reads_with_indels(); print(len(indels))
    #cov=assignment1.calculate_total_average_coverage(readlength=100);print(cov)
    #mapped=assignment1.get_number_mapped_reads();print(mapped)

    #gene_symbol=assignment1.get_gene_symbol(symbol="INS");print(gene_symbol)
    #region=assignment1.get_region_of_gene("INS");print(region)
    #exons=assignment1.get_number_of_exons("INS");print(exons)
    assignment1.print_summary()