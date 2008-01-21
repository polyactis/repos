#!/usr/bin/python

import gtk

def hello_cb(button):
	print "Hello World"
	window.destroy()

window = gtk.GtkWindow(gtk.WINDOW_TOPLEVEL) # create a top level window
window.connect("destroy", gtk.mainquit)  # quit the event loop on destruction
window.set_border_width(10)              # set padding round child widget

button = gtk.GtkButton("Hello World")
#button.connect("clicked", hello_cb)      # call hello_cb when clicked
window.add(button)                       # add button to window
button.show()                            # show button

window.show()
gtk.mainloop()                               # enter the main event loop

