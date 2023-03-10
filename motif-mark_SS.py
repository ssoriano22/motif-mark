#!/usr/bin/env python

import cairo
import math
import re

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

#Initialize pycairo canvas coordinates for display - png format
canvas_width, canvas_height = 2000, 500
surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, canvas_width, canvas_height)
#Initialize coordinates for drawing
context = cairo.Context(surface)

#Create Rectangle object
rect_1 = Rectangle(800,180,350,150)
rect_1.draw_rectangle(context)

#Create Line object 1
line_1 = Line(300,250,800,250)
line_1.draw_line(context)

#Create Line object 2
line_2 = Line(1150,250,1650,250)
line_2.draw_line(context)

#Write surface to png
surface.write_to_png("motif-mark_OoCA_SS.png")