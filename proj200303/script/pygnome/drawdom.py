#!/usr/bin/python

from gtk import *
from gnome.ui import *
import lindnaface
import GdkImlib,sys,os


class domaindraw:

	def __init__(self,fname,appwin):
		inf=open(fname,'r')
		self.appwin=appwin
		self.gl=lindnaface.infile(inf)
		self.m=lindnaface.mark(self.gl)
		self.pics=0
		self.textfont="-*-simsun-*-*-*-*-12-*-*-*-*-*-iso8859-*"
		self.picfont="-*-simsun-*-*-*-*-16-*-*-*-*-*-iso8859-*"
		self.dpainted={'Nuclear protein SET':0,
		'SET':0,
		'Nuclear protein Zn2+-binding':0,
		'Pre-SET':0,
		'YDG_SRA':0,
		'Myb DNA-binding domain':0,
		'SANT':0,
		'A+T-hook':0,
		'AT_hook':0,
		'Chromo domain':0,
		'TPR':0,
		'DUF260':0,
		'PHD':0,
		'PWWP':0,
		'zf-CXXC':0,
		'zf-MYND':0,
		'zf-C2H2':0,
		'AWS':0,
		'SET-related region':0,
		'Post-SET':0,
		'Tesmin/TSO1-like CXC domain':0,
		'FY-rich domain, N-terminal':0,
		'FY-rich domain, C-terminal':0,
		'BAH':0,
		'Whey acidic protein, core region':0,
		'Cytochrome c  heme-binding site':0,
		'ank':0,
		'Proline-rich extensin':0}
		for i in range(len(self.m)-1):
			lindnaface.sort(self.gl,self.m[i]+1,self.m[i+1])

		

	def item_event(self,widget,event=None):
		name=widget.get_data('name')
		label=self.appwin.label
		if event.type==GDK.ENTER_NOTIFY:
			label.set_text(name)
			return TRUE
		elif event.type==GDK.LEAVE_NOTIFY:
			label.set_text('')
			return TRUE
			
	
	
	def bubbledraw(self,group,i,j):
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
		'Tesmin/TSO1-like CXC domain':'RosyBrown',
		'FY-rich domain, N-terminal':'aquamarine2',
		'FY-rich domain, C-terminal':'aquamarine2',
		'BAH':'salmon',
		'WW\/Rsp5\/WWP domain':'RoyalBlue',
		'Iron hydrogenase, small subunit':'SlateGray1',
		'AP endonuclease, family 2':'red',
		'Calcium-binding EF-hand':'magenta',
		'Zinc carboxypeptidase A metalloprotease (M14)':'chocolate',
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
		'Post-SET':'rect',
		'Tesmin/TSO1-like CXC domain':'rect_r',
		'FY-rich domain, N-terminal':'tri_r',
		'FY-rich domain, C-terminal':'tri_l',
		'BAH':'elli',
		'Calcium-binding EF-hand':'dia',
		'AP endonuclease, family 2':'dia',
		'Proline-rich extensin':'tri_l'}
		
		rx1=150+int(self.gl[i][1])
		ry1=-25
		width=int(self.gl[i][2])-int(self.gl[i][1])
		height=50
	
		if width <=1:
			group.add('text',text='('+self.gl[i][2]+')',x=rx1+25,y=0,fill_color='black',font=self.textfont,anchor=ANCHOR_WEST)
			return
	
		fn=''
		for k in range(len(self.gl[i][0])):
			if self.gl[i][0][k]=='(' or self.gl[i][0][k]==')':
				pass
			else:
				fn=fn+self.gl[i][0][k]
				
		outfname='/home/gtkusr/script/bubblespng/'+fn+'_'+repr(i-self.m[j]-1)+'.png'
		print repr(i)+'\t'+repr(j)
		if os.path.isfile(outfname):
			pass
		else:
			inf=open('/tmp/bubbles.def','w')
		
			inf.write(outfname+'|')
			
			inf.write(self.gl[i][3][0:8]+'|')
			inf.write(repr(width)+'|')
			if dcolor.has_key(self.gl[i][3]):
				inf.write(dcolor[self.gl[i][3]]+'|')
			else:
				inf.write('darkslategray'+'|')
			if dshape.has_key(self.gl[i][3]):
				inf.write(dshape[self.gl[i][3]]+'\n')
			else:
				inf.write('rect'+'\n')
			
			inf.close()
			
			
			L=['run_batch.pl']
			os.spawnvp(os.P_WAIT,'./run_batch.pl',L)
		
	
		im=GdkImlib.Image(outfname)
#Next line is to work around a limitation in GdkImlib: add a reference of 'im' to an item('group')
		group.set_data('ima',im)
		ima=group.get_data('ima')
		im.render()
		widget=group.add('image',image=ima,x=rx1,y=ry1,width=width,height=height,anchor=ANCHOR_NORTH_WEST)
		#latter small emblem illustration
		pic_y=(self.pics+len(self.m)-j)*75
		if self.dpainted.has_key(self.gl[i][3]):
			if self.dpainted[self.gl[i][3]]==0:
				group.add('text',text=self.gl[i][3],x=width+25,y=pic_y+25,fill_color='black',font=self.picfont,anchor=ANCHOR_WEST)
				group.add('image',image=ima,x=0,y=pic_y,width=width,height=height,anchor=ANCHOR_NORTH_WEST)
				self.pics=self.pics+1
				print 'paint'+self.gl[i][3]
				self.dpainted[self.gl[i][3]]=1
		else:
				group.add('text',text=self.gl[i][3],x=width+25,y=pic_y+25,fill_color='black',font=self.picfont,anchor=ANCHOR_WEST)
				group.add('image',image=ima,x=0,y=pic_y,width=width,height=height,anchor=ANCHOR_NORTH_WEST)
				self.pics=self.pics+1
				print 'paint'+self.gl[i][3]
			
		widget.set_data('name','-'.join(self.gl[i]))
		widget.connect('event',self.item_event)
		
	
	
	
	def domaindraw(self,group,i):
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
		rx1=150+int(self.gl[i][1])
		ry1=-25
		rx2=150+int(self.gl[i][2])
		ry2=25
		
		# indicate the sequence length at the tail
		width=int(self.gl[i][2])-int(self.gl[i][1])
		if width <=1:
			group.add('text',text='('+self.gl[i][2]+')',x=rx1+25,y=0,fill_color='black',font=self.textfont,anchor=ANCHOR_WEST)
			return
			
		if dcolor.has_key(self.gl[i][3]):
			widget=group.add('rect',x1=rx1,y1=ry1,x2=rx2,y2=ry2,fill_color=dcolor[self.gl[i][3]],outline_color='black',width_units=1.0)
		else:
			widget=group.add('rect',x1=rx1,y1=ry1,x2=rx2,y2=ry2,fill_color='lightyellow',outline_color='blue',width_units=1.0)
	
		widget.set_data('name','-'.join(self.gl[i]))
		widget.connect('event',self.item_event)
		
		group.add('text',text=self.gl[i][3],x=(rx1+rx2)/2,y=(ry1+ry2)/2,fill_color='black',font=self.textfont,anchor=ANCHOR_CENTER,clip_width=rx2-rx1,clip_height=ry2-ry1,clip=TRUE) 


#the CORE function here

	def draw(self,canvas):
		# compute the maximum sequence length
		maxlen=0
		for i in range(len(self.gl)):
			if int(self.gl[i][2])>maxlen:
				maxlen=int(self.gl[i][2])
		
		
		gx,gy=0,0
		root=canvas.root()
		#refresh the white board
		
		root.add('rect',x1=-50,y1=-50,x2=4500,y2=2200,fill_color='white')
		#add a ruler
		group=root.add('group')
		group.add('line',points=(150,0,150+maxlen,0),width_pixels=2,fill_color='red')
		i=0
		while i*200 < maxlen-100:
			
			group.add('text',text=repr(i*200),x=150+i*200,y=0,fill_color='black',font=self.textfont,anchor=ANCHOR_SOUTH)
			i=i+1
			
		group.add('text',text=repr(maxlen),x=150+maxlen,y=0,fill_color='black',font=self.textfont,anchor=ANCHOR_SOUTH)
		
		
		j=0
		for i in range(len(self.gl)):
			if i==0:
				group=root.add('group')
				gx=gx
				gy=gy+75
				group.move(gx,gy)
				group.add('text',text=self.gl[i][0],x=0,y=0,fill_color='black',font=self.textfont,anchor=ANCHOR_WEST)
				group.add('line',points=(150,0,150+int(self.gl[i][1]),0),width_pixels=5,fill_color='slategray')
				self.bubbledraw(group,i,j)
#				self.domaindraw(group,i)				
			elif self.gl[i][0]==self.gl[i-1][0]:
				group.add('line',points=(150+int(self.gl[i-1][2]),0,150+int(self.gl[i][1]),0),width_pixels=5,fill_color='slategray')
				self.bubbledraw(group,i,j)
#				self.domaindraw(group,i)
				
			else:
				
				group=root.add('group')
				gx=gx
				gy=gy+75
				group.move(gx,gy)
				group.add('text',text=self.gl[i][0],x=0,y=0,fill_color='black',font=self.textfont,anchor=ANCHOR_WEST)
				group.add('line',points=(150,0,150+int(self.gl[i][1]),0),width_pixels=5,fill_color='slategray')
				j=j+1
				self.bubbledraw(group,i,j)
#				self.domaindraw(group,i)
	





