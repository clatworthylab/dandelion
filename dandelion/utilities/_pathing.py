#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-18 23:38:05
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-05-18 23:41:06
import dandelion
import subprocess as sp
import sys

path = dandelion.__path__[0]
command = path + "/bin"
sp.call([command] + sys.argv)