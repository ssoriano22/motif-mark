#!/usr/bin/env python

import cairo
import math
import re
import bioinfo
import argparse

class Rectangle:

    def __init__(self,x1,y1,x2,y2):
        self.start_x = x1
        self.start_y = y1
        self.end_x = x2
        self.end_y = y2
    
    def draw_rectangle(self,ctx):
        '''Draw a rectangle with pycairo. Enter start (x1,y1) and end (x2,y2) coords.'''
        ctx.rectangle(self.start_x,self.start_y,self.end_x,self.end_y) #(x1,y1,x2,y2)
        ctx.fill()
        print("Rectangle successfully drawn!")

class Line:
    
    def __init__(self,line_x1,line_y1,line_x2,line_y2):
        self.Lx1 = line_x1
        self.Ly1 = line_y1
        self.Lx2 = line_x2
        self.Ly2 = line_y2

    def draw_line(self,ctx2):
        '''Draw a line with pycairo. Enter start (x1,y1) and end (x2,y2) coords.'''
        ctx2.set_line_width(10)
        ctx2.move_to(self.Lx1,self.Ly1)  #(x1,y1)
        ctx2.line_to(self.Lx2,self.Ly2) #(x2,y2)
        ctx2.stroke()
        print("Line successfully drawn!")

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

#Argparse
def get_args():
    '''Argparse: Retrieves user-input arguments from command line.'''
    parser = argparse.ArgumentParser(description="A program to create a genetic diagram (PNG) of motif locations in a given FASTA file. Inputs: -f INPUT_FILENAME.fasta -m MOTIFS_FILENAME.txt")

    parser.add_argument("-f","--file",help="Input FASTA file",type=str)
    parser.add_argument("-m","--motifs",help="Input file with motifs",type=str)
    return parser.parse_args()

def checkMotif(seq,exonFlag):
    '''Check input sequence for specified motifs.'''
    pass

#Get argparse variables
args = get_args()
input_FASTA = args.file
input_motifs = args.motifs

#Read in motifs.txt - store motif in known_motifs_set
known_motif_set = set()
with open(input_motifs,"r") as fh_motif:
    for motif in fh_motif:
        current_motif = motif.strip()
        known_motif_set.add(current_motif)
#print(known_motif_set)


#Read in fasta and transform to onelinefasta with same name prefix
input_filename = input_FASTA.split(".fasta")[0]
output_filename = input_filename + "_temp.fasta"
bioinfo.oneline_fasta(input_FASTA,output_filename)

#Read in oneline fasta
i = 0 #Initialize counter
with open(output_filename,"r") as fh_fasta:
    while True:
        i += 1 #Increment record counter (FASTQ starts w/ 1)
        if i>1: #FOR TESTING
            break #FOR TESTING
        current_header = fh_fasta.readline().strip()
        current_genename = re.split(">",current_header)[1].split(" ")[0]
        current_geneend = re.split(">",current_header)[1].split(" ")[1].split("-")[1]
        current_genechr = re.split(">",current_header)[1].split(" ")[1].split("-")[0].split(":")[0]
        current_genestart = re.split(">",current_header)[1].split(" ")[1].split("-")[0].split(":")[1]
        current_geneseq = fh_fasta.readline().strip()
        currentGene = Gene(current_genename,current_genechr,current_genestart,current_geneend,current_geneseq)
        current_Xons_set = currentGene.get_Xons()
        #Check each sequence for each motif object, tracking intron/exon status
        for sequence in current_Xons_set:
            if sequence.isupper() == True:
                #Exon
                checkMotif(sequence,True)
            else:
                #Intron
                checkMotif(sequence,False)



#Initialize pycairo canvas coordinates for display - png format``
# canvas_width, canvas_height = 2000, 500
# surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, canvas_width, canvas_height)
# #Initialize coordinates for drawing
# context = cairo.Context(surface)

# #Create Rectangle object
# rect_1 = Rectangle(800,180,350,150)
# rect_1.draw_rectangle(context)

# #Create Line object 1
# line_1 = Line(300,250,800,250)
# line_1.draw_line(context)

# #Create Line object 2
# line_2 = Line(1150,250,1650,250)
# line_2.draw_line(context)

# #Write surface to png
# surface.write_to_png(input_filename+".png")