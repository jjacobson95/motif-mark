


import cairo

surface = cairo.SVGSurface("plot.pdf", 100, 100)
context = cairo.Context(surface)
context.set_line_width(1)
context.move_to(25,25)
context.line_to(75,25)
context.stroke()
surface.finish()
show_pdf("plot.pdf")

