import os
from PIL import Image
num=0

while num!=1:
    fileName=input("Enter image name:")
    imageName=fileName+".png"
    print(imageName)
    im=Image.open(imageName)
    im.show()
    num=input("0:rechoose, 1:continue:")
    num=int(num)

os.rename("./plot.dat","./"+fileName+"_plot.dat")
os.rename("./parameters.dat","./"+fileName+"_parameters.dat")

