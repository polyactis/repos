#!/bin/sh

#for i in $(ls /home/gtkusr/script/gp.set/gp); do
for i in `seq 1 9`; do
	echo item: gp$i
	./appdd.py /home/gtkusr/script/gp.set/gp/gp$i&
done

./appdd.py /home/gtkusr/script/gp.set/gp/gp5_5&
#./appdd.py ../SDG.xml.dom1&
#./appdd.py ../SDG.xml.dom2&
#./appdd.py ../6new.xml.dom1&


