#!/usr/bin/python

from gtk import *
from gnome.ui import *
import lindnaface
import GdkImlib,sys,os

def item_event(widget,event=None):
	name=widget.get_data('name')
	if event.type==GDK.ENTER_NOTIFY:
		label.set_text(name)
		return TRUE
	elif event.type==GDK.LEAVE_NOTIFY:
		label.set_text('')
		return TRUE
		


def bubbledraw(group,gl,i,j):
	dcolor={'Nuclear protein SET':'black',
	'SET':'black',
	'Nuclear protein Zn2+-binding':'red',
	'Pre-SET':'red',
	'YDG_SRA':'green1',
	'Myb DNA-binding domain':'aquamarine',
	'SANT':'aquamarine',
	'A+T-hook':'grey0',
	'AT_hook':'grey0',
	'Chromo domain':'pink',
	'TPR':'blue1',
	'DUF260':'magenta1',
	'PHD':'wheat',
	'PWWP':'cyan1',
	'zf-CXXC':'brown',
	'zf-MYND':'brown',
	'zf-C2H2':'brown',
	'AWS':'BlueViolet',
	'SET-related region':'turquoise1',
	'Post-SET':'turquoise1',
	'  ':'white'}
	
	dshape={'Nuclear protein SET':'pent_d',
	'SET':'pent_d',
	'Nuclear protein Zn2+-binding':'hex_h',
	'Pre-SET':'hex_h',
	'YDG_SRA':'pent_u',
	'Myb DNA-binding domain':'pent_d',
	'SANT':'pent_d',
	'A+T-hook':'rect',
	'AT_hook':'rect',
	'Chromo domain':'tri_l',
	'TPR':'rect_r',
	'DUF260':'pent_r',
	'PHD':'pent_l',
	'PWWP':'pent_r',
	'zf-CXXC':'rect',
	'zf-MYND':'rect',
	'zf-C2H2':'rect',
	'AWS':'hex_v',
	'SET-related region':'rect',
	'Post-SET':'rect'}
	
	rx1=150+int(gl[i][1])
	ry1=-25
	width=int(gl[i][2])-int(gl[i][1])
	height=50

	if width <=2:
		return

	fn=''
	for k in range(len(gl[i][0])):
		if gl[i][0][k]=='(' or gl[i][0][k]==')':
			pass
		else:
			fn=fn+gl[i][0][k]
			
	outfname='/home/gtkusr/script/bubblespng/'+fn+'_'+repr(i-m[j]-1)+'.png'
	print repr(i)+'\t'+repr(j)
	if os.path.isfile(outfname):
		pass
	else:
		inf=open('/tmp/bubbles.def','w')
	
		inf.write(outfname+'|')
		
		inf.write(gl[i][3][0:8]+'|')
		inf.write(repr(width)+'|')
		if dcolor.has_key(gl[i][3]):
			inf.write(dcolor[gl[i][3]]+'|')
		else:
			inf.write('darkslategray'+'|')
		if dcolor.has_key(gl[i][3]):
			inf.write(dshape[gl[i][3]]+'\n')
		else:
			inf.write('rect'+'\n')
		
		inf.close()
		
		
		L=['run_batch.pl']
		os.spawnvp(os.P_WAIT,'./run_batch.pl',L)
	

	im=GdkImlib.Image(outfname)
	group.set_data('ima',im)
	ima=group.get_data('ima')
	im.render()
	widget=group.add('image',image=ima,x=rx1,y=ry1,width=width,height=height,anchor=ANCHOR_NORTH_WEST)
	widget.set_data('name',gl[i][3])
	widget.connect('event',item_event)
	



def domaindraw(group,gl,i):
	dcolor={'Nuclear protein SET':'black',
	'SET':'black',
	'Nuclear protein Zn2+-binding':'red',
	'Pre-SET':'red',
	'YDG_SRA':'green1',
	'Myb DNA-binding domain':'aquamarine',
	'A+T-hook':'grey0',
	'AT_hook':'grey0',
	'Chromo domain':'pink',
	'TPR':'blue1',
	'DUF260':'magenta1',
	'PHD':'wheat',
	'PWWP':'cyan1',
	'zf-CXXC':'brown',
	'zf-MYND':'brown',
	'zf-C2H2':'brown',
	'AWS':'BlueViolet',
	'SET-related region':'turquoise1',
	'  ':'white'}
	rx1=150+int(gl[i][1])
	ry1=-25
	rx2=150+int(gl[i][2])
	ry2=25
	if dcolor.has_key(gl[i][3]):
		group.add('rect',x1=rx1,y1=ry1,x2=rx2,y2=ry2,fill_color=dcolor[gl[i][3]],outline_color='black',width_units=1.0)
	else:
		group.add('rect',x1=rx1,y1=ry1,x2=rx2,y2=ry2,fill_color='lightyellow',outline_color='azure1',width_units=1.0)

	group.add('text',text=gl[i][3],x=(rx1+rx2)/2,y=(ry1+ry2)/2,fill_color='black',font='-*-*-*-*-*--*-*-*-*-*-60-*-*',anchor=ANCHOR_CENTER,clip_width=rx2-rx1,clip_height=ry2-ry1,clip=TRUE)


def draw(gl,canvas):
	gx,gy=0,0
	root=canvas.root()
	j=0
	for i in range(len(gl)):
		if i==0:
			group=root.add('group')
			group.add('text',text=gl[i][0],x=0,y=0,fill_color='black',font='-*-*-*-*-*--*-*-*-*-*-60-*-*',anchor=ANCHOR_NORTH_WEST)
			group.add('line',points=(150,0,150+int(gl[i][1]),0),width_pixels=5,fill_color='slategray')
			bubbledraw(group,gl,i,j)
			
		elif gl[i][0]==gl[i-1][0]:
			group.add('line',points=(150+int(gl[i-1][2]),0,150+int(gl[i][1]),0),width_pixels=5,fill_color='slategray')
			bubbledraw(group,gl,i,j)
			
		else:
			group=root.add('group')
			gx=gx
			gy=gy+100
			group.move(gx,gy)
			group.add('text',text=gl[i][0],x=0,y=0,fill_color='black',font='-*-*-*-*-*--*-*-*-*-*-60-*-*',anchor=ANCHOR_NORTH_WEST)
			group.add('line',points=(150,0,150+int(gl[i][1]),0),width_pixels=5,fill_color='slategray')
			j=j+1
			bubbledraw(group,gl,i,j)






if len(sys.argv)>1:
	infile=sys.argv[1]
else:
	infile=raw_input("please enter the Pfam file:\n")

inf=open(infile,'r')

gl=lindnaface.infile(inf)

m=lindnaface.mark(gl)

for i in range(len(m)-1):
	lindnaface.sort(gl,m[i]+1,m[i+1])

		

win=GtkWindow()
sw=GtkScrolledWindow()
sw.set_policy(POLICY_AUTOMATIC,POLICY_AUTOMATIC)

win.connect('destroy',mainquit)
win.set_title('Domain Drawing Program')

vbox=GtkVBox()
win.add(vbox)
vbox.show()

label=GtkLabel()
vbox.pack_start(label,expand=FALSE)

canvas=GnomeCanvas()
canvas.set_scroll_region(-50,-100,4500,1300)
canvas.set_pixels_per_unit(0.5)
canvas.root().add('rect',x1=-50,y1=-100,x2=4500,y2=1300,fill_color='white')

sw.add(canvas)
vbox.add(sw)
canvas.show()
draw(gl,canvas)
win.set_default_size(600,550)
win.show_all()
mainloop()
