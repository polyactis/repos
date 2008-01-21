#!/usr/bin/env python2.2

# This is a sample implementation of an editor.

from gtk import *
import GtkExtra
from gnome.ui import *
import GdkImlib,sys,os

#My modules
import lindnaface,drawdom



class AppWindow(GtkWindow):
	def __init__(self, quit_cb=None):
		GtkWindow.__init__(self, WINDOW_TOPLEVEL)
		self.set_usize(600, 500)
		self.connect("delete_event", self.file_exit)
		
		self.quit_cb = quit_cb
		self.imgfname=''
		self.remember_x=0
		self.remember_y=0
		self.new_x=0
		self.new_y=0
		self.ratio=1
		
		self.table = GtkTable(3,1)
		self.add(self.table)
		self.table.show()
		hbox=GtkHBox()
		self.table.attach(hbox, 0,1, 0,1, xoptions=FILL,
				  yoptions=FILL)

		hbox.show()
		self.menubar = self.create_menu()
		hbox.add(self.menubar)
		self.menubar.show()

		self.label=GtkLabel()
		self.table.attach(self.label,0,1,1,2,yoptions=SHRINK)


		self.sw=GtkScrolledWindow()
		self.sw.set_policy(POLICY_AUTOMATIC,POLICY_AUTOMATIC)
		self.table.attach(self.sw,0,1,2,3)

		self.canvas=GnomeCanvas()
#		self.canvas.connect('event',self.canvas_event)
		self.canvas.set_scroll_region(-50,-50,4500,2200)
		self.canvas.set_pixels_per_unit(self.ratio)
#		self.canvas.root().add('rect',x1=-50,y1=-50,x2=4500,y2=1300,fill_color='white')

		self.sw.add(self.canvas)

	def canvas_event(self,widget,event=None):
		if event.type==GDK.BUTTON_PRESS:
		#	if event.state & GDK.BUTTON1_MASK:
				self.remember_x=event.x*2-50
				self.remember_y=event.y*2-50
		elif event.type==GDK.MOTION_NOTIFY:
			if event.state & GDK.BUTTON1_MASK:
				self.new_x=event.x*2-50
				self.new_y=event.y*2-50
		
		elif event.type==GDK.BUTTON_RELEASE:
			widget.root().add('line',points=(self.remember_x,self.remember_y,self.new_x,self.remember_y),width_pixels=1)
			widget.root().add('line',points=(self.remember_x,self.remember_y,self.remember_x,self.new_y),width_pixels=1)
			widget.root().add('line',points=(self.remember_x,self.new_y,self.new_x,self.new_y),width_pixels=1)
			widget.root().add('line',points=(self.new_x,self.remember_y,self.new_x,self.new_y),width_pixels=1)
			self.remember_x=0
			self.remember_y=0
			self.new_x=0
			self.new_y=0


	def load_file(self, fname):
#		try:

			dd=drawdom.domaindraw(fname,self)
			dd.draw(self.canvas)
			self.set_title(os.path.abspath(fname))
			self.imgfname=os.path.abspath(fname)+'.png'
			self.fname = fname

#		except:
#			GtkExtra.message_box('DomainDraw', "Can't open " + fname,("OK",))

	def create_menu(self):
		mf = GtkExtra.MenuFactory()

		mf.add_entries([
			('File/Open...',    '<control>O', self.file_open),
			('File/<separator>',None,         None),
			('File/New Window','<control>N',self.new_window),
			('File/<separator>',None,         None),
			('File/Exit',       '<control>Q', self.file_exit),
			('Tool/Zoom In','=',self.zoom_in),
			('Tool/Zoom Out','-',self.zoom_out),
			('Save Image','<control>S',self.save_image),
			('Help/About...',   None,         self.help_about)
		])
		# activate key bindings ...
		self.add_accel_group(mf.accelerator)
		self.mf = mf
		return mf
	
	def zoom_in(self,mi=None):
		self.ratio=self.ratio+0.1
		self.canvas.set_pixels_per_unit(self.ratio)
		self.canvas.update_now()
	
	def zoom_out(self,mi=None):
		self.ratio=self.ratio-0.1
		self.canvas.set_pixels_per_unit(self.ratio)
		self.canvas.update_now()

	def fs_close(self,b):
		self.fs.hide()
		self.fs.destroy()
		
	def save_image(self,mi=None):
		try:
			window=self.canvas.get_window()		
			self.image=GdkImlib.create_image_from_drawable(window,None,-50,-100,4500,1800)
			self.image.save_image(self.imgfname)
		except:
			GtkExtra.message_box('No Image available',("OK",))
		
	def new_window(self,mi=None):
		w=AppWindow(quit_cb=mainquit)
		w.show_all()
	
		
	def file_open(self, mi=None):
		fname = GtkExtra.file_open_box(modal=FALSE)
		if not fname: return
		self.load_file(fname)

	def file_exit(self, mi=None, event=None):
		self.hide()
		
		self.destroy()
		if self.quit_cb: self.quit_cb(self)

	def help_about(self, mi):
		GtkExtra.message_box("Domain Draw Window", "Copyright (C) 2003  " +
				     "Yu Huang\n" +
				     "This program is covered by the GPL.",
				     ("OK",))




	
if __name__ == '__main__':
	import sys

	
	quit_cb = mainquit

	w = AppWindow(quit_cb=quit_cb)
	
	if len(sys.argv)>1:
		w.load_file(sys.argv[1])
			
	w.show_all()
	w.set_usize(0,0)
	mainloop()


