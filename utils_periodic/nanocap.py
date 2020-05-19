#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from forestUtils import Forest
import pickle
from atomUtils import Atom
# In[ ]:


"body class"

class Body(object):
    def __init__(self):
        
        #x,y,z,r represent vectors containing x,y,z locations of the points
        # rvector contains radius of the cylinder on which x,y points are located
        
        self._x, self._y, self._z, self._r = None, None, None, None
        self._rfname, self._sl = None, None # file name and starting line from which
        self._N = 0 # number of atoms in the body
        
        self._id, self._t, self._q = None, None, None
    
    # properties or attributes of the class Body 
    #(easy to access these attributes this way)
    @property
    def N(self):
        return len(self.x)
    
    @property
    def sl(self):
        return self._sl
    @property
    def x(self):
        return self._x
    @property
    def y(self):
        return self._y
    @property
    def z(self):
        return self._z
    @property
    def r(self):
        return self._r
    @property
    def rfname(self):
        return self._rfname
    @property
    def id(self):
        return self._id
    @property
    def t(self):
        return self._t
    @property
    def q(self):
        return self._q
    
    
    

    
    
    # setter functions for attributes of the class Body
    # easy to set values this way
    @sl.setter
    def sl(self,data):
        self._sl = data
        return
    @x.setter
    def x(self,data):
        self._x = data
        return
    @y.setter
    def y(self,data):
        self._y = data
        return
    @z.setter
    def z(self,data):
        self._z = data
        return
    @r.setter
    def r(self,data):
        self._r = data
        return
    @rfname.setter
    def rfname(self,data):
        self._rfname = data
        return
    

    
    
    # read information from .xyz file and store the values in respective attributes
    def read_xyz(self):
        
        rfile = open(self.rfname,'r')
        lcount = 0
        #initialization
        xv,yv,zv,rv = [],[],[],[]
        
        #read line by line
        for line in rfile:
            if lcount >= self.sl:
                tmp = line.split()
                #convert from string to float
                x,y,z = float(tmp[1]), float(tmp[2]), float(tmp[3])
                #compute r
                r = np.sqrt(x**2+y**2)
                
                #append the x,y,z,r vectors
                xv.append(x), yv.append(y) 
                zv.append(z), rv.append(r)
            
            lcount += 1
        
        # update the attributes of the class Body
        self.x, self.y, self.z, self.r = xv, yv, zv, rv
        
        return None


    def read_data(self):

        # read lammps file

        rfile = open(self.rfname,'r')
        lcount = 0

        idv, tv, qv = [],[],[]
        xv, yv, zv,rv = [],[],[],[]

        for line in rfile:
            if lcount >= sl:
                tmp = linesplit()

                sno, t, q = float(tmp[0]), float(tmp[1]), float(tmp[2])
                x, y, z = float(tmp[3]), float(tmp[4]), float(tmp[5])

                r = np.sqrt(x**2 + y**2)

                idv.append(sno), tv.append(t). qv.append(q)
                xv.append(x), yv.append(y), zv.append(z)
                rv.append(r)

            lcount += 1

        self.id, self.t, self.q = idv, tv, qv
        self.x, self.y, self.z = xv, yv, zv
        self.r = rv 

        return None

        


# In[ ]:


"cap class"

class Cap(object):
    
    #this class defines attributes of a Cap for a nanotube
    
    def __init__(self):
        
        # x,y,z,r are vectors containing coordinates and the radius of the sphere
        # on which a particular atom is located
        self._x, self._y, self._z, self._r = None, None, None, None
        
        self._rfname, self._sl = None, None # attributes of the file from which info is read
        self._N = 0 # number of atoms in the Cap

        self._id, self._t, self._q = None, None, None

        return
    
    
    # properties or attributes of the class Cap
    # easy to access this way
    @property
    def N(self):
        return len(self.x)
    @property
    def sl(self):
        return self._sl
    @property
    def x(self):
        return self._x
    @property
    def y(self):
        return self._y
    @property
    def z(self):
        return self._z
    @property
    def r(self):
        return self._r
    @property
    def rfname(self):
        return self._rfname
    @property
    def id(self):
        return self._id
    @property
    def t(self):
        return self._t
    @property
    def q(self):
        return self._q
    
    
    
    
    
    # setter functions for all the attributes of the class Cap
    # easy to set values this way
    @sl.setter
    def sl(self,data):
        self._sl = data
        return
    @x.setter
    def x(self,data):
        self._x = data
        return
    @y.setter
    def y(self,data):
        self._y = data
        return
    @z.setter
    def z(self,data):
        self._z = data
        return
    @r.setter
    def r(self,data):
        self._r = data
        return
    @rfname.setter
    def rfname(self,data):
        self._rfname = data
        return
    
    
    
    # function to read information from a file and store them as 
    # attributes of the class
    
    def read_xyz(self):
        
        rfile = open(self.rfname,'r')
        lcount = 0
        xv,yv,zv,rv = [],[],[],[]
        for line in rfile:
            if lcount >= self.sl:
                tmp = line.split()
                x,y,z = float(tmp[1]), float(tmp[2]), float(tmp[3])
                r = np.sqrt(x**2 + y**2 + z**2)
                
                xv.append(x), yv.append(y), zv.append(z),rv.append(r)
            
            lcount += 1
        
        
        self.x, self.y, self.z, self.r = xv,yv,zv,rv
        
        return None
                


# In[ ]:


"class cnt"

class Cnt(object):
    

    def __init__(self):
        
        # attributes of this class are
        # body
        # primary cap
        # secondary cap
        # each being a class of its own
        
        self._body, self._primaryCap, self._secondaryCap = Body(), Cap(), Cap()
        self._N = None

        return
    
    @property
    def body(self):
        return self._body
    
    @property
    def primaryCap(self):
        return self._primaryCap
    
    @property
    def secondaryCap(self):
        return self._secondaryCap
    @property
    def N(self):
        return self.body.N + self.primaryCap.N + self.secondaryCap.N
    

    @body.setter
    def body(self,data):
        self._body = data
        return
    
    @primaryCap.setter
    def primaryCap(self,data):
        self._primaryCap = data
        return
    
    @secondaryCap.setter
    def secondaryCap(self,data):
        self._secondaryCap = data
        return 
    
    def read_body(self,fname,sl):
        
        # read information about body of a cnt from file fname
        # store information in the respective attribute
        
        self.body.rfname = fname
        self.body.sl = sl 
        self.body.read_xyz()
        
        return
    
    def read_cap(self,fname,sl, flag=''):
        
        # read information about a cap of a cnt from file fname
        # if flag == primary, its primary or bottom cap
        # if flag == secondary, its secondary or top cap
        
        if flag == 'primary':

            self.primaryCap.rfname = fname
            self.primaryCap.sl = sl
            self.primaryCap.read_xyz()

        if flag == 'secondary':

            self.secondaryCap.rfname = fname
            self.secondaryCap.sl = sl
            self.secondaryCap.read_xyz()

        
        return
    
    def scale_body_h(self,cnt2):
        
        # scale height of the body of the cnt to equal that of cnt2
        
        # computing Hmax and Hmin of the cnt2
        Hmax, Hmin = max(cnt2.body.z), min(cnt2.body.z)
        
        # computing hmax and hmin of the current cnt
        hmax, hmin = max(self.body.z), min(self.body.z)
        
        # computing shift on current cnt and ratio
        alpha, r = hmin - Hmin, Hmax/hmax
        
        # initializing an empty list
        zvec = []
        
        # iterating through all the atoms of current cnt body and shifting and scaling them
        for i in range(self.body.N):
            tmp = (self.body.z[i]-alpha)*r
            zvec.append(tmp)
        
        return zvec
    
    def scale_body_r(self,cnt2):
        
        # scale radius of the body of cnt to equal that of cnt2
        
        # computing max radius of cnt2 and current cnt
        R, r = max(cnt2.body.r), max(self.body.r)
        
        #initializing empty vectors
        xvec,yvec,rvec = [],[],[]
        
        # computing the ratio to which the current locations of the atoms should be scaled to
        ratio = R/r
        
        # iterating through all the atoms of current cnt body and scaling 
        # x,y values and recomputing new r
        
        for i in range(self.body.N):
            xnew,ynew = self.body.x[i]*ratio, self.body.y[i]*ratio
            xvec.append(xnew), yvec.append(ynew)
            
            rnew = np.sqrt(xnew**2 + ynew**2)
            rvec.append(rnew)
            
        return xvec,yvec, rvec
    
    
    def scale_body(self,cnt2):
        
        # scale body of the cnt to match that of cnt2
        
        #initializing new instance of Body
        bodyNew = Body()
        
        # computing new x,y, locations scaled according to radius of cnt2
        bodyNew.x, bodyNew.y, bodyNew.r = self.scale_body_r(cnt2)
        
        # computing new z locations scaled according height of cnt2
        bodyNew.z = self.scale_body_h(cnt2)
        
        return bodyNew
        
    
    def scale_cap(self,R,zshift):


        # computing Hmax and Hmin of the cnt2
        # Hmin = max(cnt2.body.z), min(cnt2.body.z)
        
        # # computing hmax and hmin of the current cnt
        # hmin = max(self.body.z), min(self.body.z)
        
        # # computing shift on current cnt and ratio
        # alpha = hmin - Hmin
        
        newPrimaryCap = Cap()
        newSecondaryCap = Cap()
        
        rPrimary, rSecondary = max(self.primaryCap.r), max(self.secondaryCap.r)
        
        ratioPrimary, ratioSecondary = R/rPrimary, R/rSecondary
        
        xvec, yvec, zvec,zvecSecondary, rvec = [],[],[],[],[]
        
        for i in range(self.primaryCap.N):


            
            xnew = self.primaryCap.x[i]*ratioPrimary
            ynew = self.primaryCap.y[i]*ratioPrimary
            znew = self.primaryCap.z[i]*ratioPrimary
            znewSecondary = -znew + zshift

            rnew = np.sqrt(xnew**2 + ynew**2 + znew**2)
                
            xvec.append(xnew), yvec.append(ynew), zvec.append(znew), rvec.append(rnew)
            zvecSecondary.append(znewSecondary)
        
        newPrimaryCap.x, newPrimaryCap.y, newPrimaryCap.z, newPrimaryCap.r = xvec, yvec, zvec, rvec
        
        
        newSecondaryCap.x, newSecondaryCap.y, newSecondaryCap.z, newSecondaryCap.r = xvec, yvec, zvecSecondary, rvec
        
        return newPrimaryCap, newSecondaryCap
        
        
         
    
    def scale_cnt(self,cnt2):
        
        newBody = self.scale_body(cnt2)
        
        RnewCap = max(newBody.r)
        zshift = max(newBody.z)

        newPrimaryCap, newSecondaryCap = self.scale_cap(RnewCap,zshift)
        
        scaledCnt = Cnt()
        scaledCnt.body = newBody
        scaledCnt.primaryCap = newPrimaryCap
        scaledCnt.secondaryCap = newSecondaryCap
        
        # return newBody, newPrimaryCap, newSecondaryCap
        return scaledCnt


    def write_xyz(self,dir = '', fname = '', body=True, primaryCap = True, secondaryCap = True):

        
        if body:

            wfile = open(dir + 'carbon_lattice_Nanotube_scaled.xyz','w')
            wfile.write(str(self.body.N) + '\n\n')

            for i in range(self.body.N):

                tmp = ['C', str(self.body.x[i]), str(self.body.y[i]), str(self.body.z[i])]
                wfile.write(' '.join(tmp) + '\n')

            wfile.close()

        if primaryCap:

            wfile = open(dir + 'carbon_lattice_Cap Primary_scaled.xyz','w')
            wfile.write(str(self.primaryCap.N) + '\n\n')

            for i in range(self.primaryCap.N):

                tmp = ['C', str(self.primaryCap.x[i]), str(self.primaryCap.y[i]), str(self.primaryCap.z[i])]
                wfile.write(' '.join(tmp) + '\n')

            wfile.close()

        if secondaryCap:

            wfile = open(dir + 'carbon_lattice_Cap Secondary_scaled.xyz','w')
            wfile.write(str(self.secondaryCap.N) + '\n\n')

            for i in range(self.secondaryCap.N):

                tmp = ['C', str(self.secondaryCap.x[i]), str(self.secondaryCap.y[i]), str(self.secondaryCap.z[i])]
                wfile.write(' '.join(tmp) + '\n')


            wfile.close()



        return None

    def write_data(self,dirName='',fName='',flag='cap',q=0.0):


        wfile = open(dirName+fName,'w')

        if flag == 'cap':
            N = self.body.N + self.secondaryCap.N 
        elif flag == 'body':
            N = self.body.N


        wfile.write('# LAMMPS data file written by OVITO\n')
        wfile.write(str(N) + ' atoms\n')
        wfile.write('1 atom types\n\n')
        wfile.write('Atoms # charge\n\n')

        lcount = 1

        for i in range(self.body.N):
            tmp = [lcount,1,q,self.body.x[i],self.body.y[i],self.body.z[i]]

            wfile.write(' '.join(list(map(str,tmp))) + '\n')

            lcount += 1

        if flag == 'cap':

            for i in range(self.secondaryCap.N):

                tmp = [lcount, 1, q, self.secondaryCap.x[i], self.secondaryCap.y[i], self.secondaryCap.z[i]]

                wfile.write(' '.join(list(map(str,tmp))) + '\n')

                lcount += 1

        wfile.close()

        return None

    def write_data_grid(self,dirName='',fName='',flag='cap',q=0.0,G=[2,2,25,25]):

        xGrid, yGrid = G[0], G[1]
        xLambda, yLambda = G[2], G[3]

        wfile = open(dirName+fName,'w')

        lcount = 0

        if flag == 'cap':
            N = self.body.N + self.secondaryCap.N
        elif flag == 'body':
            N = self.body.N

        N = N*xGrid*yGrid

        wfile.write('# LAMMPS data file written by OVITO\n')
        wfile.write(str(N) + ' atoms\n')
        wfile.write('1 atom types\n\n')
        wfile.write('Atoms # charge\n\n')

        ncount = 1

        for xind in range(xGrid):

            for yind in range(yGrid):

                xshift = xind*xLambda
                yshift = yind*yLambda

                for i in range(self.body.N):

                    x = self.body.x[i] + xshift
                    y = self.body.y[i] + yshift
                    z = self.body.z[i]


                    tmp = [ncount,1,q,x,y,z]

                    wfile.write(' '.join(list(map(str,tmp))) + '\n')

                    ncount += 1

                if flag =='cap':

                    for i in range(self.secondaryCap.N):

                        x = self.secondaryCap.x[i] + xshift
                        y = self.secondaryCap.y[i] + yshift
                        z = self.secondaryCap.z[i]

                        tmp = [ncount, 1, q, x, y, z]

                        wfile.write(' '.join(list(map(str,tmp))) + '\n')

                        ncount += 1

        wfile.close()  
        return None
        
    
    def save(self,fname='',q = 0.0,G = [1,1,0,0],bflag=True,scapflag=True,pcapflag=False):

        nx, ny = G[0],G[1]
        lx, ly = G[2],G[3]

        fst = Forest()
        fst.nx, fst.ny = nx, ny
        fst.lx, fst.ly = lx, ly

        unit = Cnt()
        unit.body = self.body
        unit.primaryCap = self.primaryCap
        unit.secondaryCap = self.secondaryCap

        fst.unit = unit
        fst.bodyr = max(self.body.r)
        fst.bodyh = max(self.body.z)
        fst.h = max(self.secondaryCap.z)


        N = bflag * self.body.N + pcapflag * self.primaryCap.N + scapflag * self.secondaryCap.N

        N = N * nx * ny

        ncount = 1

        atoms = []

        for xind in range(nx):

            for yind in range(ny):

                xshift = xind*lx
                yshift = yind*ly

                for i in range(self.body.N):

                    x = self.body.x[i] + xshift
                    y = self.body.y[i] + yshift
                    z = self.body.z[i] 

                    newAtom = Atom()
                    newAtom.x = x
                    newAtom.y = y
                    newAtom.z = z
                    newAtom.q = q
                    newAtom.id = ncount
                    newAtom.type = 1
                    atoms.append(newAtom)

                    ncount += 1

                if scapflag :

                    for i in range(self.secondaryCap.N):

                        x = self.secondaryCap.x[i] + xshift
                        y = self.secondaryCap.y[i] + yshift
                        z = self.secondaryCap.z[i]

                        newAtom = Atom()
                        newAtom.x = x
                        newAtom.y = y 
                        newAtom.z = z 
                        newAtom.q = q 
                        newAtom.id = ncount
                        newAtom.type = 1
                        atoms.append(newAtom)

                        ncount += 1

        fst.atoms = atoms
        
        
        

        pickle.dump(fst,open(fname,'wb'))


        return None

    

