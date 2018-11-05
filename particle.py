from PIL import Image
import os
import errno

fileName=input("Enter file name.:")
imgName=fileName
fileName=fileName+".png"
im=Image.open(fileName)#image file name for initial condition

rgb_im=im.convert("RGB")

size=rgb_im.size

rgb_im2=rgb_im.transpose(Image.FLIP_TOP_BOTTOM)#to make lower left origin
iBP=0
iFLP=0
iOBP=0
f=open("wall.txt","w")
f2=open("numbers.h","w")
for x in range(size[0]):
    for y in range(size[1]):
        r,g,b=rgb_im2.getpixel((x,y))
        s="%d %d %d %d %d \n"%(x,y,r,g,b)
        if (r,b)==(127,127):
            f.write(s)
            iBP=iBP+1

f.close()            
f=open("fluid.txt","w")
for x in range(size[0]):
    for y in range(size[1]):
        r,g,b=rgb_im2.getpixel((x,y))
        if(g==97):
            print(r,g,b)
        s="%d %d %d %d %d \n"%(x,y,r,g,b)
        if (r,g,b)==(127,127,255):
            f.write(s)
            iFLP=iFLP+1
       

f.close()            
f=open("obstacle.txt","w")
for x in range(size[0]):
    for y in range(size[1]):
        r,g,b=rgb_im2.getpixel((x,y))
        s="%d %d %d %d %d \n"%(x,y,r,g,b)
        if (r,b)==(255,127):
            f.write(s)
            iOBP=iOBP+1
#SPH.c中にOPBが分母になっている箇所があり、OBP=0ではエラーになってしまう。
#そのためOBP=0の時、計算の範囲外にOBPを一つ配置しておく。
if(iOBP==0):
    s="%d %d %d %d %d \n"%(-10,-10,255,100,127)
    f.write(s)
    iOBP=iOBP+1
f.close()

s2="//%s\n"%(fileName)
f2.write(s2)
s2="#define FLP %d \n"%(iFLP)
f2.write(s2)
s2="#define BP %d \n"%(iBP)
f2.write(s2)
s2="#define OBP %d \n"%(iOBP)
f2.write(s2)
s2="#define N %d \n"%(iFLP+iBP+iOBP)
f2.write(s2)
s2="#define XSIZE %d\n"%(size[0])
f2.write(s2)
s2="#define YSIZE %d\n"%(size[1])
f2.write(s2)

f2.close()

s2='Source_%s'%(imgName)
print(s2)
try:
    os.mkdir(s2)
except OSError as e:
    if e.errno==errno.EEXIST:
        print("Directory already exists")
    else:
        raise




