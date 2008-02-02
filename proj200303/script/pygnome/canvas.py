#! /usr/bin/env python

from gtk import *
import gnomecanvas
from gnome.ui import *
#import GdkImlib

win = Window()
#button=GtkButton('EMBEDDED')
win.connect('destroy', main_quit)
win.set_title('Canvas test')

canvas = gnomecanvas.Canvas(aa=True)
#canvas.set_usize(300, 300)
canvas.set_scroll_region(0,0,300,300)
win.add(canvas)

#board=[]
#canvas.set_data('board',board)
root=canvas.root()
group=root.add(gnomecanvas.CanvasGroup())
#board.append(group)
group.add(gnomecanvas.CanvasLine,points=(0,0, 100,100),width_pixels=10,fill_color='slategray')

#group=GnomeCanvasGroup()
#group=root.add('group')
#gp=root.get_data('group')
#group=canvas.root().add('group')
#board.append(group)
group.move(50,50)
group.add(gnomecanvas.CanvasPolygon, points=(10,10, 90,10, 90,90, 10,90, 0,50),
		  width_units=5, fill_color='slategray')
#group.add('rect',x1=100,y1=100,x2=200,y2=150,fill_color='lightyellow',outline_color='black',width_units=1.0)

w = canvas.root().add(gnomecanvas.CanvasText,text='SET',x=150,y=125,fill_color='black',font='-adobe-helvetica-bold-r-normal--*-240-75-75-p-*-*-1',anchor=ANCHOR_CENTER)
w.set(weight=800)
canvas.root().add(gnomecanvas.CanvasText, text='SET', x=200, y=175, fill_color='red', anchor=ANCHOR_E)


#file='/tmp/bubbles/B_lectin.png'
#im=GdkImlib.Image(file)
#im.render()
#group.add('image',image=im,x=20,y=200,width=180,height=50,anchor=ANCHOR_NORTH_WEST)

#group.add('widget',widget=button,x=5,y=5)
win.show_all()


main()

