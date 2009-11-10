#!/usr/bin/env python
import GADA
ins = GADA.GADA()
segments_ls = ins.run([1,1,1,1,0.99,0.99,1,1,0.1,0.1,0.1,0.15], 0.2, 4, 2)
print segments_ls
