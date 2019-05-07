from ovito.io import import_file
from ovito.io import *
from ovito.vis import *
from ovito.modifiers import *
from math import *
import sip
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import numpy as np
from numpy import linalg as LA
from ovito.vis import Viewport
import sys

filename = sys.argv[1]


pipeline = import_file(filename, columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z"], multiple_frames = True)
pipeline.add_to_scene()

print(filename)


type_list = pipeline.source.data.particles.particle_types.types

radius = 0.56
for index in range(len(type_list)):
    type_list[index].radius = radius

patch_types = [1,2,7,8]
dark_purple = (0.635294118, 0.0, 1.0)
bright_purple = (0.894117647, 0.698039216, 0.905882353)
for index in range(len(type_list)):
    if index in patch_types:
        type_list[index].color = dark_purple
    else:
        type_list[index].color = bright_purple

pipeline.source.data.cell.vis.render_cell = False

vp = Viewport(type = Viewport.Type.Front)
vp.zoom_all()
#vp.camera_pos = (0,250,0)
vp.camera_dir=(200,200,200)








vp.render_image(filename = './fibrils.png',size=(1000,1000),frame=1,background=(0.,0.,0.),alpha=False,renderer=TachyonRenderer(ambient_occlusion=True,ambient_occlusion_brightness=0.7,ambient_occlusion_samples=16,antialiasing=True,antialiasing_samples=16,direct_light=True,direct_light_intensity=2,shadows=True))
