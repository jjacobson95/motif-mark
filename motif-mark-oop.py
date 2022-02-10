#!usr/bin/env/python3.9

#download cairo from conda forge. use new envrionment for this.

import cairo
#import math

#0,0 is in top left corner of plot.
surface = cairo.SVGSurface("plot.svg", 100, 100)
context = cairo.Context(surface)
context.set_line_width(1)
context.move_to(25,25)
context.line_to(75,25)
context.stroke()
context.finish()
show_svg("plot.svg")
