#!/usr/bin/env python

import cairo
import math
import re
import bioinfo
import argparse

class Exon:

    def __init__(self,exon_seq,x1,y1,x2,y2):
        self.eseq = exon_seq
        self.start_x = x1
        self.start_y = y1
        self.end_x = x2
        self.end_y = y2
    
    def draw_exon(self,ctx):
        '''Draw a exon (thicker line) with pycairo. Enter start (x1,y1) and end (x2,y2) coords.'''
        #ctx.rectangle(self.start_x,self.start_y,self.end_x,self.end_y) #(x1,y1,x2,y2)
        ctx.set_line_width(100)
        ctx.move_to(self.start_x,self.start_y)  #(x1,y1)
        ctx.line_to(self.end_x,self.end_y) #(x2,y2)
        ctx.stroke()
        print("Exon successfully drawn!")

class Intron:
    
    def __init__(self,intron_seq,line_x1,line_y1,line_x2,line_y2):
        self.iseq = intron_seq
        self.Lx1 = line_x1
        self.Ly1 = line_y1
        self.Lx2 = line_x2
        self.Ly2 = line_y2

    def draw_intron(self,ctx2):
        '''Draw a intron (line) with pycairo. Enter start (x1,y1) and end (x2,y2) coords.'''
        ctx2.set_line_width(10)
        ctx2.move_to(self.Lx1,self.Ly1)  #(x1,y1)
        ctx2.line_to(self.Lx2,self.Ly2) #(x2,y2)
        ctx2.stroke()
        print("Intron successfully drawn!")

class Gene:

    def __init__(self,gene_Name,gene_Chr,gene_Start,gene_End,gene_Seq):
        self.name = gene_Name
        self.chr = gene_Chr
        self.startLoc = gene_Start
        self.endLoc = gene_End
        self.seq = gene_Seq

    def get_Xons(self):
        '''Parse input fasta sequence into exon (uppercase seq), and upstream/downstream introns (lowercase seq). Returns a tuple(intron_upstream,exon,intron_downstream).'''
        intron_up = re.split("[A-Z]+",self.seq)[0]
        intron_down = re.split("[A-Z]+",self.seq)[1]
        exon = re.split("[a-z]+",self.seq)[1]
        Xons = (intron_up, exon, intron_down)
        return Xons
    
class Motif:

    def __init__(self, motif_seq):
        self.mseq = motif_seq
    
    def searchMotif(self,cseq,exonFlag):
        '''Check input sequence for specified motifs.'''
        mlength = len(self.mseq)
        if "Y" in self.mseq:
            #Ambiguous bases in motif
            pass
        else:
            #No ambiguous bases in motif
            pass

#Argparse
def get_args():
    '''Argparse: Retrieves user-input arguments from command line.'''
    parser = argparse.ArgumentParser(description="A program to create a genetic diagram (PNG) of motif locations in a given FASTA file. Inputs: -f INPUT_FILENAME.fasta -m MOTIFS_FILENAME.txt")

    parser.add_argument("-f","--file",help="Input FASTA file",type=str)
    parser.add_argument("-m","--motifs",help="Input file with motifs",type=str)
    return parser.parse_args()

#Get argparse variables
args = get_args()
input_FASTA = args.file
input_motifs = args.motifs

#Read in motifs.txt - store motif in known_motifs_set
known_motif_set = set()
with open(input_motifs,"r") as fh_motif:
    for motif in fh_motif:
        current_motif = Motif(motif.strip())
        known_motif_set.add(current_motif)
#print(known_motif_set)


#Read in fasta and transform to onelinefasta with same name prefix
input_filename = input_FASTA.split(".fasta")[0]
output_filename = input_filename + "_temp.fasta"
bioinfo.oneline_fasta(input_FASTA,output_filename)

#Initialize pycairo canvas coordinates for display - png format``
canvas_width, canvas_height = 2000,1000
surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, canvas_width, canvas_height)
#Initialize coordinates for drawing
context = cairo.Context(surface)

#Read in oneline fasta
i = 0 #Initialize counter - for testing
n = 10 #Gene counter - to adjust x coords for drawing per gene
with open(output_filename,"r") as fh_fasta:
    while True:
        i += 1 #Increment record counter (starts w/ 1)
        if i>1: #FOR TESTING - one gene
            break #FOR TESTING
        #Get current gene header and split into name, chr, start/end loci
        current_header = fh_fasta.readline().strip()
        #EOF case:
        if current_header == "":
            break
        current_genename = re.split(">",current_header)[1].split(" ")[0]
        current_geneend = re.split(">",current_header)[1].split(" ")[1].split("-")[1]
        current_genechr = re.split(">",current_header)[1].split(" ")[1].split("-")[0].split(":")[0]
        current_genestart = re.split(">",current_header)[1].split(" ")[1].split("-")[0].split(":")[1]
        #Get current gene sequence + create Gene object
        current_geneseq = fh_fasta.readline().strip()
        currentGene = Gene(current_genename,current_genechr,current_genestart,current_geneend,current_geneseq)
        #Break current sequence into introns and exon. Format: tuple(up_intron,exon,down_intron)
        current_Xons = currentGene.get_Xons()
        #Create Intron/Exon objects for each of current gene (2 Intron, 1 Exon)
        currentUpIntron = Intron(current_Xons[0],n,500,(n+len(current_Xons[0])),500)
        currentExon = Exon(current_Xons[1],(n+len(current_Xons[0])),500,(n+len(current_Xons[0]+current_Xons[1])),500)
        currentDownIntron = Intron(current_Xons[2],(n+len(current_Xons[0]+current_Xons[1])),500,(n+len(current_Xons[0]+current_Xons[1]+current_Xons[2])),500)
        #Draw introns and exon objects to context/pycairo surface
        currentUpIntron.draw_intron(context)
        currentExon.draw_exon(context)
        currentDownIntron.draw_intron(context)

        #Check each sequence for each motif object, tracking intron/exon status
        for mot in known_motif_set:
            for sequence in current_Xons:
                if sequence.isupper() == True:
                    #Exon
                    mot.searchMotif(sequence,True)
                else:
                    #Intron
                    mot.searchMotif(sequence,False)
        #Increment n according to total x length drawn for this gene
        n += len(current_Xons[0]+current_Xons[1]+current_Xons[2])+10 #remove +10 to connect each gene

#Write surface to png - name using same prefix as input FASTA file
surface.write_to_png(input_filename+".png")
print("Motif Mark PNG generated successfully! Image Name:",input_filename+".png")