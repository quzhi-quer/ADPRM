import pandas as pd
import numpy as np
import struct
import sys
from os      import system
import matplotlib.pyplot as plt 
import py3Dmol
from IPython.display import display,HTML
#import mdtraj
import re
#required files: mol.csv contains element mass/charge informations
def readpdb(fi='target.pdb',libfile='mol.csv'):
    pdb=open(fi,'r')
    index=[]
    name=[]
    resname=[]
    chain=[]
    resid=[]
    x=[]
    y=[]
    z=[]
    charge=[]
    element=[]
    mass=[]
    atomdict={}
    lib=pd.read_csv(libfile)
    lib.element=lib.element.fillna('NA')
    for i in lib.index:
        atomdict[lib.element[i]]=[lib.mass[i]]
    lines=pdb.readlines()
    print(atomdict)
    for line in lines:
        if line[0:4] == 'ATOM':
            index.append(int(line[6:11]))
            n=line[12:16]
            name.append(n.strip())
            resname.append(line[17:20])
            resid.append(int(line[22:26]))
            x.append(float(line[30:38]))
            y.append(float(line[38:46]))
            z.append(float(line[46:54]))
            if len(line)<77 or line[76:78].strip()=='':
                element.append(n.strip()[0])
            else:
                element.append(line[76:78].strip())
            mass.append(atomdict[element[-1]][0])
#            charge.append(atomdict[element[-1]][1])
    data=pd.DataFrame({
        'series':index,
        'name':name,
        'resname':resname,
        'resid':resid,
        'x':x,
        'y':y,
        'z':z,
        'mass':mass,
#        'charge':charge,
        'element':element
    })
    data[['series', 'resid']] = data[['series', 'resid']].astype(int)
    pdb.close()
    return data
def writepdb(data,outfile='out.pdb'):
    pdb=open(outfile,'w')
    out=[]
    for i in data.index:
        out.append("ATOM  {0:>5} {1:<4} {2:<3} {3}{4:>4}    {5:>8.3f}{6:>8.3f}{7:>8.3f}                      {8:>2}\n"
              .format(int(data.series[i]),data.name[i],data.resname[i],'X',int(data.resid[i]),data.x[i],data.y[i],data.z[i],data.element[i]))
    pdb.writelines(out)
    pdb.close()
    return
def readtrj(inpfile='target.trj',ref=[],snap=(0,2500,1),name='traj'):
    data=open(inpfile,'rb')
    data.read(81)
    pos=np.zeros((len(ref.index),int((snap[1]-snap[0])/snap[2]),3),dtype=np.float16)
    for s in range(snap[0],snap[1],snap[2]):
        count=0
#        print('read frame ',s)
        for i in range(0,len(ref.index)):
            for x in range(0,3):
                pos[i][s][x]=float(data.read(8))
                count+=1
                if count>=10:
                    data.read(1)
                    count=0
        data.read(26)
        ref[name+'x'+str(s)]=pos[:,s,0]
        ref[name+'y'+str(s)]=pos[:,s,1]
        ref[name+'z'+str(s)]=pos[:,s,2]
    data.close()
    return ref
def writetrj(data,outfile='out.trj',trajname='traj',snap=(0,2500,1),title='Amber trjectory file write by quzhi'):
    fo=open(outfile,'wb')
    fo.write(("{0:<80}".format(title)).encode())
    fo.write(b'\x0A')
    for s in range(snap[0],snap[1],snap[2]):
        x=trajname+'x'+str(s)
        y=trajname+'y'+str(s)
        z=trajname+'z'+str(s)
        count=0
        for i in data.index:
            fo.write("{0:>8.3f}".format(data[x][i]).encode())
            count+=1
            if count>=10:
                count=0
                fo.write(b'\x0A')
            fo.write("{0:>8.3f}".format(data[y][i]).encode())
            count+=1
            if count>=10:
                count=0
                fo.write(b'\x0A')
            fo.write("{0:>8.3f}".format(data[z][i]).encode())
            count+=1
            if count>=10:
                count=0
                fo.write(b'\x0A')
        fo.write(b'\x0A')
        fo.write("{0:>8.3f}".format(9999).encode())
        fo.write("{0:>8.3f}".format(9999).encode())
        fo.write("{0:>8.3f}".format(9999).encode())
        fo.write(b'\x0A')
    fo.close()
    return
def writeparm(data,parmfile='out.parm'):
#this parm file do not contain parameters!!!!!
#only use for curve analysis
    fp=open(parmfile,'w')
    fp.writelines('''%VERSION  VERSION_STAMP = V0001.000  DATE = 10/09/21 20:21:36                  
%FLAG TITLE                                                                     
%FORMAT(20a4)                                                                   
default_name                                                                    
%FLAG POINTERS                                                                  
%FORMAT(10I8)
''')
    fp.writelines("{0:>8d}{1:>8d}{2:>8d}{3:>8d}{4:>8d}{5:>8d}{6:>8d}{7:>8d}{8:>8d}{9:>8d}\n".\
    format(len(data.index),len(data.name.value_counts()),0,0,0,0,0,0,0,0))
#  NATOM    : total number of atoms 
#  NTYPES   : total number of distinct atom types
    fp.writelines("{0:>8d}{1:>8d}{2:>8d}{3:>8d}{4:>8d}{5:>8d}{6:>8d}{7:>8d}{8:>8d}{9:>8d}\n".\
    format(0,len(data.resid.value_counts()),0,0,0,0,0,0,0,0))
#resid now only support amber format pdb
#
    fp.writelines('''%FLAG ATOM_NAME                                                                 
%FORMAT(20a4)                                                                                                                               
''')
    line=''
    count=0
    for i in data.index:
        line=line+"{0:<4}".format(data.name[i])
        count+=1
        if count>=20:
            count=0
            fp.writelines(line+'\n')
            line=''
    fp.writelines(line+'\n')
    fp.writelines('''%FLAG RESIDUE_LABEL                                                             
%FORMAT(20a4)                                                                   
''')
    line=''
    count=0
    temp=''
    respointer=[]
    for i in data.index:
        if data.resid[i]!=temp:
            respointer.append(i+1)
            temp=data.resid[i]            
            line=line+"{0:<4}".format(data.resname[i])
            count+=1
            if count>=20:
                count=0
                fp.writelines(line+'\n')
                line=''
    fp.writelines(line+'\n')
    fp.writelines('''%FLAG RESIDUE_POINTER                                                           
%FORMAT(10I8)                                                                   
''') 
    line=''
    count=0
    for i in respointer:
        line=line+"{0:>8d}".format(i)
        count+=1
        if count>=10:
            count=0
            fp.writelines(line+'\n')
            line=''
    fp.writelines(line+'\n')
    fp.close()
    return
def select(data,selname='sel',mode=1,relation=1,sel={'index':[],'series':[],'name':[],'resname':[],'resid':[],'chain':[],'element':[]}):
#mode 1/0 represent select/deselect
#relation 1/0 represent and/or
    if relation==0:
        data[selname]=0
        for name in sel:
            if sel[name]:
                data[selname][data[name] in sel[name]]=mode
    else:
        data[selname]=1
        for name in sel:
            if sel[name]:
                for i in data.index:
                    if data[name][i] not in sel[name]:
                        data[selname][i]=1-mode
    return
def addcom(data,selname='sel',trajname='traj',weight='mass'):
    select=data[data[selname]==1]
    totalmass=select[weight].sum()
    com={'series':0,'resname':selname,'resid':len(data.resid.value_counts())+1,
         'mass':totalmass,'element':'X','name':selname,
         'x':select[weight].dot(select.x)/totalmass,
         'y':select[weight].dot(select.y)/totalmass,
         'z':select[weight].dot(select.z)/totalmass}
    for s in data:
        if s.startswith(trajname):            
            com[s]=select[weight].dot(select[s])/totalmass
    return data.append(com,ignore_index=True)
def genlib(data,selname='',flib='/home/pangpang/curve+/standard',pre='X',charge=1,mapfile='temp.csv'):
    system('mv '+flib+'_i.lib backup')
    fo=open(flib+'_i.lib','w')
    fm=open()
    count=0
    output=data['index','name','resid','resname',selname].copy
    for i in output.index:
        if output[selname][i]:
            fo.writelines(pre+"{:<3d}".format(count)+' '+str(charge)+'\n')
            fm.writelines(str(int(output.index[i]))+','+output.name[i]+','
                          +str(int(output.resid[i]))+','+output.resname[i]+','+pre+"{:<3d}".format(count)+',\n')
            count+=1
    fo.close()
    fm.close()
    return output

def curves(data,sel,deffnm='curve_dna',inp='',outdir='/home/pangpang',trajname='traj',
           snap=(0,800,1),libori='/home/pangpang/curve+/standard',pre='X',IFwritetrj=True):
#
    datap=data[['name','resname','resid']].copy()
    if inp=='':
        fi=open(deffnm+'.inp','w')
        fi.writelines('''export inp='''+deffnm+'''.trj
export parm='''+deffnm+'''.prmtop 
export out='''+deffnm+'''
/home/pangpang/curve+/Cur+ <<EOF
 &inp file=$inp, ftop=$parm, lis=$out, ions=.t., frames=.t.,
 lib='''+deffnm+''', &end
2 1 -1 0 0
1:20
40:21
EOF''')
        fi.close()
    print('sh command: cp '+libori+'_b.lib '+deffnm+'_b.lib')
    system('cp '+libori+'_b.lib '+deffnm+'_b.lib')
    print('sh command: cp '+libori+'_s.lib '+deffnm+'_s.lib') 
    system('cp '+libori+'_s.lib '+deffnm+'_s.lib')    
    fl=open(deffnm+'_i.lib','w')            
#standard_i.lib file in curves is now a map file write
    count=1
    for i in data[data[sel]==1].index:
        fl.writelines(pre+str(count)+'   '+str(i)+' \n')
        datap.loc[i,'name']=pre+str(count)
        count+=1
    fl.close()
    writeparm(datap,deffnm+'.prmtop')
    writetrj(data,outfile=deffnm+'.trj',trajname=trajname,snap=snap,title='taskname:'+deffnm+' for curves analysis')
    system('sh '+deffnm+'.inp')  
    readcdi(deffnm+'.cdi',ref=data).to_csv(deffnm+'.csv')
    system('mv '+deffnm+'* '+outdir)
    return 

def readcdi(fi,ref=''):
    cdi=open(fi,'r',encoding='utf-8')
    lines=cdi.readlines()
    index=lines[2].split()
    ionlist=lines[3].split()
    frames=int((len(lines)-5)/2)
    ions=len(index)
    table={}
    ionlist={}
    circle={}
    if type(ref)==pd.DataFrame:    
        c=1
        for i in index:
            ionlist[c]=ref['resname'][int(i)]#+str(ref['resid'][int(i)])+ref['name'][int(i)]
            circle[ionlist[c]]=0
            c+=1
    else:
        c=1
        for i in index:
            ionlist[c]=index
            c+=1
    for i in ionlist:
        table[ionlist[i]+'pos']=np.zeros(frames)
        table[ionlist[i]+'rad']=np.zeros(frames)
        table[ionlist[i]+'ang']=np.zeros(frames)
    f=-1
    ln=4
    print(frames)
    while ln<len(lines)-1:
        l=0
        f+=1
        ln+=1        
        ni=lines[ln]
        if ni.startswith('      0'):
            continue
        else:
            ni=int(ni)
            ln+=1
        line=lines[ln]
        for i in range(ni):
            ion=ionlist[int(line[l+26:l+28])]
            table[ion+'pos'][f]=float(line[l:l+7])/1000
            table[ion+'rad'][f]=float(line[l+7:l+14])/1000
            table[ion+'ang'][f]=float(line[l+14:l+21])*0.057324841
            l+=28
    curvesdata=pd.DataFrame(table)
    return curvesdata
def dataref(data,maxfix=10,county=1000,maxrot=120):
#step I 
#replace zeros to nan
#del collums with data less than county
#then del rows from beginning&end until met rows with whole data
    for i in data:
        if i.endswith('pos'):
            data[i][data[i]==0]=np.nan
    for i in data:
        if i.endswith('pos') and data.count()[i]<county:
            data=data.drop([i],axis=1)
            data=data.drop([i[:-3]+'rad'],axis=1)
            data=data.drop([i[:-3]+'ang'],axis=1)
            print('del collums:',i[:-3])
    outrange=[0,data.shape[0]-1]
    while True:
        switch=0
        for i in data:
            if i.endswith('pos') and np.isnan(data[i][outrange[0]]):
                outrange[0]+=1
                switch=1
                break
        if switch==0:
            break
    while True:
        switch=0
        for i in data:
            if i.endswith('pos') and np.isnan(data[i][outrange[1]]):
                outrange[1]-=1
                switch=1
                break
        if switch==0:
            break
    print('select data:',outrange)
    data=data[outrange[0]:outrange[1]+1].copy()    
#step II fix data with zero position values
    data.reset_index(drop=True)
    lmax=len(data)
    added=0
    for i in data:
        if i.endswith('pos')==False:
            continue
        f=0
        while f<lmax:
            if np.isnan(data[i][f])==False:
                up=[data[i][f],data[i[:-3]+'rad'][f],data[i[:-3]+'ang'][f]]
                down=(0,0,0)
                f+=1
            else:
                count=0
                while(np.isnan(data[i][f])):    
                    count+=1
                    f+=1
                if count>maxfix:
                    raise ValueError("too many missing data in row",i,f)
                down=[data[i][f],data[i[:-3]+'rad'][f],data[i[:-3]+'ang'][f]]
                if down[2]-up[2]>180:
                    up[2]+=360
                elif down[2]-up[2]<-180:
                    up[2]-=360
                s=[(down[0]-up[0])/(count+1),(down[1]-up[1])/(count+1),(down[2]-up[2])/(count+1)]
                for c in range(1,count+1):
                    data[i][f-c]=down[0]-c*s[0]
                    data[i[:-3]+'rad'][f-c]=down[1]-c*s[1]
                    data[i[:-3]+'ang'][f-c]=down[2]-c*s[2]
                    print('add data ',i[:-3],f-c,data[i][f-c],data[i[:-3]+'rad'][f-c],data[i[:-3]+'ang'][f-c])
                    added+=1
    print('total data added:',added)
#step III fix data with rotation more than 180 degree between frames
    warn=0
    for i in data:
        if i.endswith('ang')==False:
            continue
        encircle=0
        for f in range(1,lmax):
            da=data[i][f]-data[i][f-1]-360*encircle
            if da%360>maxrot and da%360<(360-maxrot):
                warn+=1
                data[i][f]=data[i][f-1]
                print('WARNING',warn,': ',i,' rotated ',da,' degree at frame ',f)
                continue
            elif da>maxrot:
#                print(i[:-3],' encircle +1 with DNA at frame ',f)
                encircle+=1
            elif da<-maxrot:
#                print(i[:-3],' encircle -1 with DNA at frame ',f)
                encircle-=1
            data[i][f]-=360*encircle
    return data
def dataunity(data):
    posstd=0
    radstd=0
    angstd=0
    for i in data:
        if i.endswith('pos'):
            posstd+=data[i].std()
        elif i.endswith('rad'):
            radstd+=data[i].std()
        elif i.endswith('ang'):
            angstd+=data[i].std()
    for i in data:
        if i.endswith('pos'):
            data[i]/=(posstd*3/len(data.columns))
        elif i.endswith('rad'):
            data[i]/=(radstd*3/len(data.columns))
        elif i.endswith('ang'):
            data[i]/=(angstd*3/len(data.columns))
    return data
def data_average(data,window=100):
    ld=len(data)-window
    adata=data[:1-window].copy(deep=True)
    for i in adata:
        temp=(data[i][window:].values-data[i][:-window].values)/window
        adata[i][0]=adata[i][:window].mean()
        for f in range(0,ld):
            adata[i][f+1]=adata[i][f]+temp[f]
#        print(i)
    return adata
def data_diff(data,diff=1):
    ddata=data[:-diff].copy(deep=True)
    for i in ddata:
        ddata[i]=data[i][:-diff].values-data[i][diff:].values
    return ddata
def genheatmat(datax,datay,hx=10,hy=10):
    if len(datax)!=len(datay):
        print('two lines do not match!')    
    histx=(datax.max()-datax.min())/hx
    histy=(datay.max()-datay.min())/hy
    minx=datax.min()
    miny=datay.min()
    lx=len(datax)
    heatmap=np.zeros([hx+1,hy+1])
    for i in range(lx):
        heatmap[int((datax[i]-minx)//histx)][int((datay[i]-miny)//histy)]+=1
    return heatmap
