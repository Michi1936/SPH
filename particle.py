from PIL import Image

fileName=input("Enter file name.:")
fileName=fileName+".png"
im=Image.open(fileName)#image file name for initial condition

rgb_im=im.convert("RGB")

size=rgb_im.size

rgb_im2=rgb_im.transpose(Image.FLIP_TOP_BOTTOM)#to make lower left origin
iBP=0
iFLP=0
iOBP=1
f=open("wall.txt","w")
f2=open("numbers.h","w")
for x in range(size[0]):
    for y in range(size[1]):
        r,g,b=rgb_im2.getpixel((x,y))
        s="%d %d %d %d %d \n"%(x,y,r,g,b)
        if (r,g,b)==(127,127,127):
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
        if (r,g,b)==(127,97,255):
            f.write(s)
            iFLP=iFLP+1

f.close()            
f=open("obstacle.txt","w")
for x in range(size[0]):
    for y in range(size[1]):
        r,g,b=rgb_im2.getpixel((x,y))
        s="%d %d %d %d %d \n"%(x,y,r,g,b)
        if (r,g,b)==(255,0,0):
            f.write(s)
            iOBP=iOBP+1
            
s2="#define FLP %d \n"%(iFLP)
f2.write(s2)
s2="#define BP %d \n"%(iBP)
f2.write(s2)
s2="#define OBP %d \n"%(iOBP)
f2.write(s2)
s2="#define N %d \n"%(iFLP+iBP+iOBP)
f2.write(s2)
f.close()
f2.close()
