# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 12:40:33 2020

@author: lixiangk
"""

def astra2transport(infile='ast.0480.001',outfile='ast.0480.002'):
    ''' generate a new file: x/cm, xp/mrad, y/cm, yp/mrad, z/cm, Pc/GeV 
        return beam momentum Pc/GeV and normalized sigma matrix''' 
    data=np.loadtxt(infile)
    data[0,2]=0
    #data[1:,2]=data[1:,2]+data[0,2]
    data[1:,5]*=1
    data[1:,5]=data[0,5]+data[1:,5]

    data[:,0]=data[:,0]*100*1 # x: m->cm
    data[:,1]=data[:,1]*100*1 # y: m->cm
    data[:,2]=data[:,2]*100 # z: m->cm

    data[:,6]=np.sqrt(data[:,3]**2+data[:,4]**2+data[:,5]**2)/1e9 # Pc=P*c: eV->GeV
    data[:,3]=data[:,3]/data[:,5]*1e3*1 # xp: mrad
    data[:,4]=data[:,4]/data[:,5]*1e3*1 # yp: mrad
    data[:,5]=data[:,6] # pz->Pc

    data[:,6]=data[:,1] # 6th column: y
    data[:,1]=data[:,3] # 2nd column: xp

    data[:,3]=data[:,4] # 4th column: yp
    data[:,4]=data[:,2] # 5th column: z
    data[:,2]=data[:,6] # 3rd column: y

    #np.savetxt('ast.0480.002',data[:,0:6],fmt='%12.4E')
    Pc=np.mean(data[:,5])
    data[:,5]=data[:,5]/Pc*100.

    sigma=np.cov(data[:,0:6].T)
    for i in np.arange(6):
        for j in np.arange(i):
            sigma[i,j]=sigma[i,j]/np.sqrt(sigma[i,i]*sigma[j,j])
            sigma[j,i]=sigma[i,j]
            #print str.format('%12.4f' % sigma[i,j]),
        #print str.format('%12.4f' % np.sqrt(sigma[i,i]))
    for i in np.arange(6):
        sigma[i,i]=np.sqrt(sigma[i,i])
    return [Pc, sigma]
#astra2transport()

class Transport:
    def __init__(self):
        self.Pc=0
        self.sigma=np.zeros((6,6))
        self.title='/Optics of xxx facility/\n'+'0\n'+'13. 48. /ANGL/ ;\n\n'
        self.sentinel='\nSENTINEL\n/*PLOT*/\n-1\n\n'
        self.filename='for001.dat'
        self.fitting=''
        self.content=''
    
    def beam_from_astra(self,infile):
        self.Pc, self.sigma=astra2transport(infile)

    def add_title(self,title=''):
        if title=='':
            title=self.title
        self.content=self.content+title
    def add_sentinel(self,sentinel=''):
        if sentinel=='':
            sentinel=self.sentinel
        self.content=self.content+sentinel
    
    def add_text(self,text='\n'):
        self.content=self.content+text

    def add_beam(self,label='BEAM'):
        oo='1.0 '
        for i in np.arange(6):
            oo=oo+' '+str.format('%.6f' % self.sigma[i,i])
        oo=oo+' '+str.format('%.6f' % self.Pc)
        oo=oo+' /'+label+'/ ;'+'\n';
        self.content=self.content+oo
    def add_beam_corr(self,label='CORR'):
        oo='12. '
        for i in np.arange(1,6):
            for j in np.arange(i):
                oo=oo+' '+str.format('%.6f' % self.sigma[i,j])
        oo=oo+' /'+label+'/ ;'+'\n\n'
        self.content=self.content+oo
    def add_beam_type(self,i,value):
        #oo='16. '+`i`+' '+`value`+' ;\n'
        oo = str.format('16. %d %g ;\n' % (i, value))
        self.content=self.content+oo

    def add_rotation(self,beta,label=''):
        #oo='2.0'+' '+`beta`+' /'+label+'/ ;'+'\n'
        #oo='2.0'+' '+`beta`+' ;\n'
        oo = str.format('2.0 %g ;\n' % beta)
        self.content=self.content+oo
        
    def add_rotation_beam(self,alpha,label=''):
        #oo='20.0'+' '+`alpha`+' /'+label+'/ ;'+'\n'
        #oo='20.0'+' '+`alpha`+' ;\n'
        oo = str.format('20.0 %g ;\n' % alpha)
        self.content=self.content+oo
        
    def add_drift(self,length,label):
        oo='3.0'+' '+str.format('%.6f' % length)+' /'+label+'/ ;'+'\n'
        self.content=self.content+oo
    def add_bend(self,length,alpha,strength,label):
        oo='4.0'+' '+str.format('%.6f %.6f %.6f /%s/ ;\n' % (length,alpha,strength, label))
        self.content=self.content+oo
    def add_quad(self,length,strength,bore,label):
        oo='5.0'+' '+str.format('%.6f %.6f %g /%s/ ;\n' % (length,strength, bore, label))
        self.content=self.content+oo
    
    def add_special(self, digit, value):
        oo='16.'+' '+str.format('%.0f %.6f ;\n' % (digit, value))
        self.content=self.content+oo
        
    def add_update_R1(self):
        oo='6. 0. 1. ;\n'
        self.content=self.content+oo
    def add_update_R2(self):
        oo='6. 0. 2. ;\n'
        self.content=self.content+oo
    
    def add_fit(self,n,i,j,value,accuracy=0.001,label='F001'):
        #oo='-10.'+`n`+' '+`i`+' '+`j`+' '+str.format('%.6f %.6f' % (value,accuracy))+' /'+label+'/ ;'+'\n'
        oo=str.format('-10.%d %d %d %.6f %.6f' % (n, i, j, value, accuracy))+' /'+label+'/ ;'+'\n'
        #self.fitting=self.fitting+'10.'+`n`+' /'+label+'/ ;'+'\n'
        self.content=self.content+oo
    def add_fit_emittance(self,n,value,accuracy=0.001,label='XEMI'):
        oo='24. 1. 1. 1. ; \n\
24. 2. 2. 2. ; \n\
24. 12 1. 3. ; \n\
24. 12. 1. 4. ; \n\
24. 100. 1. 5. ; \n\
25. 3. 4. 3. 6. ; \n\
25. 5. 6. 2. 3. ; \n\
25. 3. 1. 5. 4. ; \n\
25. 1. 2. 3. 3. ; \n\
25. 3. 4. 3. 1. ; \n\
-10.'+str.format('%g' % n)+' 9. 1. '+str.format('%.6f %.6f' % (value,accuracy))+' /'+label+'/ ;'+'\n'
        #self.fitting=self.fitting+'10.'+`n`+' /'+label+'/ ;'+'\n'
        self.content=self.content+oo
    def add_fitting(self,n, label):
        #oo='10.'+`n`+' /'+label+'/ ;'+'\n'
        oo=str.format('10.%d /%s/ ;\n' % (n, label))
        self.content=self.content+oo
        
    def add_vary_code(self,code,flag,label):
        #oo=`code`+'.'+flag+' /'+label+'/ ;'+'\n'
        oo = str.format('%d.%s /%s/ ;\n' % (code, flag, label))
        self.content=self.content+oo
        
    def output(self,filename=''):
        if filename=='':
            filename=self.filename
        #f = file(filename,'w') # open for 'w'riting
        #f.write(self.content)  # write text to file
        #f.close()   # close the file
        
        with open(filename, 'w') as f_handle:
            f_handle.write(self.content+'\n')
        
#%% for PITZ Flash RT
from universal import *
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2020\BioRad\Case2/ARC'
os.chdir(workdir)

tran = Transport()
tran.beam_from_astra('ast.2550.007')

gamma = np.sqrt(1+(tran.Pc*1e3/g_mec2)**2)
xmax, ymax = 0.5, 0.5
xmatch, ymatch = 0.065, 0.065
zmatch = 0.235
xemi = 3/gamma


Lquad, Rquad = 0.0675, 2.15

Q01, Q02 = 0.5, -0.1
D01, D02, D03=5.6635-4.8-Lquad/2.,0.25-Lquad,6.05-5.6635-Lquad/2-0.25+(58.9+5.5)/1000.
#print D01,D02,D03

B1 = tran.Pc*1e9/600e-3/g_c
beta = 0  #17.355

Q1, Q2 = 1, -0.5

L, rho, theta = 2.230, 0.3, np.pi/3
gap = 5.0 # cm
beta1, beta2 = 30, 30

ub = (L/2-rho*(1-np.cos(theta)))/np.sin(theta)-Lquad*2;  print(ub)
#ub = 0.5
#ub = 0.384615
#ub = 0.8


D1, D2, D3 = 0.3, 0.2, 0.1
D1 = ub-D3-D2+0
#D3 = 0.5-D1-D2


Q11, Q12 = 0.0, -0.0
D11, D12, D13 = 0.2, 0.1-Lquad, 0.1

tran.add_title()
tran.add_beam()
tran.add_beam_corr()
tran.add_beam_type(31., 0.)
tran.add_beam_type(22., 50000)
tran.add_text()

#tran.add_drift(D01, 'D01')
#tran.add_quad(Lquad, Q01, Rquad, 'Q01')
#tran.add_drift(D02, 'D02')
#tran.add_quad(Lquad, Q02, Rquad, 'Q02')
#tran.add_drift(D03, 'D03')
#tran.add_text()

tran.add_drift(0.5, 'D')

# tran.add_text()
# tran.add_fit(2, 1, 1, xmax, label = 'FIT1')
# tran.add_fit(2, 3, 3, ymax, label = 'FIT2')
# tran.add_text()

tran.add_update_R1()
tran.add_text()

Ldip = rho*theta
tran.add_rotation(beta1, 'R11')
tran.add_special(5, gap/2.0)
tran.add_bend(Ldip, theta/np.pi*180, 0, 'B1') # no transverse gradient
#tran.add_bend(0.62831, 60, B1, 'B1')
tran.add_rotation(beta2, 'R12')
tran.add_text()

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D1, 'D1')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_quad(Lquad, Q1, Rquad, 'Q1')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D2, 'D2')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_quad(Lquad, Q2, Rquad, 'Q2')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D3, 'D3')

tran.add_fit(0, -1, 6, 0, label = 'R160')
tran.add_fit(0, -1, 1, 0, label = 'R110')
tran.add_fit(0, -2, 2, 0, label = 'R220')
tran.add_fit(0, -3, 3, 0, label = 'R330')
tran.add_fit(0, -4, 4, 0, label = 'R440')

tran.add_text()

tran.add_drift(D3, 'D3')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_quad(Lquad, Q2, Rquad, 'Q2')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D2, 'D2')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_quad(Lquad, Q1, Rquad, 'Q1')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D1, 'D1')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_rotation_beam(180)
tran.add_rotation(beta2, 'R21')
#tran.add_bend(0.62831, 60, B1, 'B1')
tran.add_bend(Ldip, theta/np.pi*180, 0, 'B1')
tran.add_rotation(beta1, 'R22')
tran.add_rotation_beam(-180)

#tran.add_text()
#tran.add_fit(2, 1, 1, xmax, label = 'FIT1')
#tran.add_fit(2, 3, 3, ymax, label = 'FIT2')
#tran.add_text()

tran.add_fit(0, -1, 6, 0, label = 'R160')
tran.add_fit(0, -2, 6, 0, label = 'R260')

#tran.add_fit(0, 1, 1, xmatch, label = 'FIT3')
#tran.add_fit(0, 3, 3, ymatch, label = 'FIT4')
#tran.add_fit(2, 5, 5, zmatch, label = 'FIT5')
#tran.add_text()

tran.add_fit(0, 1, 1, xmatch, label = 'F110')
tran.add_fit(0, 3, 3, ymatch, label = 'F330')
tran.add_fit(0, 5, 5, zmatch, label = 'F550')
tran.add_text()

for _ in np.arange(1):
    tran.add_drift(0.5, 'D')
    tran.add_fit(0, -1, 6, 0, label = 'FI16')
    tran.add_fit(0, -2, 6, 0, label = 'R260')
    tran.add_text()
    
tran.add_fit(0, 1, 1, xmatch, label = 'F110')
tran.add_fit(0, 3, 3, ymatch, label = 'F330')
tran.add_fit(0, 5, 5, zmatch, label = 'F550')
tran.add_text()


# tran.add_quad(Lquad, Q11, Rquad, 'Q11')
# tran.add_drift(D12, 'D12')
# #tran.add_fit(0, -1, 6, 0, label = 'FI16')
# tran.add_fit(0, 1, 1, xmatch, label = 'FIT3')
# tran.add_fit(0, 3, 3, ymatch, label = 'FIT4')
# tran.add_fit(2, 5, 5, zmatch, label = 'FIT5')
# tran.add_text()

# tran.add_drift(D12, 'D12')
# tran.add_quad(Lquad, Q12, Rquad, 'Q12')
# #tran.add_fit(0, -1, 6, 0, label = 'FI16')
# tran.add_fit(0, 1, 1, xmatch, label = 'FIT3')
# tran.add_fit(0, 3, 3, ymatch, label = 'FIT4')
# tran.add_fit(2, 5, 5, zmatch, label = 'FIT5')
# tran.add_text()

# tran.add_drift(0.1, 'D13')
# #tran.add_fit(0, -1, 6, 0, label = 'FI16')
# tran.add_fit(0, 1, 1, xmatch, label = 'FIT3')
# tran.add_fit(0, 3, 3, ymatch, label = 'FIT4')
# tran.add_fit(2, 5, 5, zmatch, label = 'FIT5')
# tran.add_text()


tran.add_fit_emittance(2, xemi)
tran.add_sentinel()

#tran.add_vary_code(5,'0F','Q01')
#tran.add_vary_code(5,'0G','Q02')
tran.add_vary_code(5, '0K', 'Q1')
tran.add_vary_code(5, '0L', 'Q2')

tran.add_vary_code(3, '-B', 'D1')
tran.add_vary_code(3, 'B', 'D2')

# tran.add_vary_code(5,'0M','Q11')
# tran.add_vary_code(5,'0N','Q12')


tran.add_text()

tran.add_fitting(2, 'F112')
tran.add_fitting(2, 'F332')
tran.add_fitting(0, 'F110')
tran.add_fitting(0, 'F330')
tran.add_fitting(0, 'F550')

tran.add_fitting(0, 'R160')
tran.add_fitting(0, 'R260')

# tran.add_fitting(0, 'R110')
# tran.add_fitting(0, 'R220')
# tran.add_fitting(0, 'R330')
# tran.add_fitting(0, 'R440')

tran.add_fitting(2,'XEMI')
tran.add_text()

tran.add_sentinel('\nSENTINEL\nSENTINEL\n')

tran.output('C:\\TRANS\\for001.dat')
#os.system('cp for001.dat C:\\TRANS\\for001.dat')


#%% For FLASH RT Chengdu
from universal import *
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2020\FRT\_B1-274Gs-B2-180Gs'
os.chdir(workdir)

tran = Transport()
tran.beam_from_astra('ast.0750.012')

gamma = np.sqrt(1+(tran.Pc*1e3/g_mec2)**2)


xmax, ymax = 0.5, 0.5
xmatch, ymatch = 0.05, 0.05
zmatch = 0.05
xemi = 10/gamma

xmax, ymax = 0.5, 0.5
xmatch, ymatch = 0.1, 0.1
zmatch = 0.075
xemi = 10/gamma


Lquad, Rquad = 0.0675, 2.15

Q01, Q02 = 0, 0
D01, D02, D03 = 0.15-Lquad/2., 0.20106-Lquad, 0.15-Lquad/2. 
print(D01, D02, D03)

#
rho, theta = 0.2, np.pi/4
gap = 5.0 # cm

B1 = tran.Pc*1e9/rho/g_c
beta = 0  #17

#
Q11, Q12, Q13 = 1, -0.5, 1
#Q11, Q12, Q13 = 1, 0., 1
#Q11, Q12, Q13 = 0, -1.0, 0
D11, D12, D13, D14 = 0.3-Lquad/2.0, 0.2-Lquad, 0.2-Lquad, 0.3-Lquad/2.0
D11, D12, D13, D14 = 0.364435, 0.3, 0.3, 0.364435

#
B2 = tran.Pc*1e9/rho/g_c
beta = 0

#
Q21, Q22, Q23 = 1, -0.25, 1
D21, D22, D23, D24 = 0.3-Lquad/2.0, 0.2-Lquad, 0.2-Lquad, 0.3-Lquad/2.0
D21, D22, D23, D24 = 0.3, 0.2, 0.2, 0.3


#
B3 = tran.Pc*1e9/rho/g_c
beta = 0


tran.add_title()
tran.add_beam()
tran.add_beam_corr()
tran.add_beam_type(31., 0.)
tran.add_beam_type(22., 50000)
tran.add_text()

tran.add_drift(D01, 'D01')
tran.add_quad(Lquad, Q01, Rquad, 'Q01')
tran.add_drift(D02, 'D02')
tran.add_quad(Lquad, Q02, Rquad, 'Q02')
tran.add_drift(D03, 'D03')
tran.add_text()

# tran.add_text()
# tran.add_fit(2, 1, 1, xmax, label = 'FIT1')
# tran.add_fit(2, 3, 3, ymax, label = 'FIT2')
# tran.add_text()

tran.add_update_R1()
tran.add_text()

theta = np.pi/4-6.504661e-2
Ldip = rho*theta
beta1, beta2 = 0, 0

tran.add_rotation(beta1, 'R11')
tran.add_special(5, gap/2.0)
tran.add_bend(Ldip, theta/np.pi*180, 0, 'B1') # no transverse gradient
#tran.add_bend(0.62831, 60, B1, 'B1')
tran.add_rotation(beta2, 'R12')
tran.add_text()

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D11, 'D11')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_quad(Lquad, Q11, Rquad, 'Q11')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D12, 'D12')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_quad(Lquad/2, Q12, Rquad, 'Q12')
#tran.add_fit(0, -1, 6, 0, label = 'R160')
#tran.add_fit(0, -2, 6, 0, label = 'R260')
tran.add_quad(Lquad/2, Q12, Rquad, 'Q12')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D12, 'D12')

tran.add_quad(Lquad, Q13, Rquad, 'Q13')

tran.add_drift(D14, 'D11')

# tran.add_text()
# tran.add_fit(2, 1, 1, xmax, label = 'F112')
# tran.add_fit(2, 3, 3, ymax, label = 'F332')
# tran.add_text()

# tran.add_fit(0, -1, 6, 0, label = 'R160')
# tran.add_fit(0, -1, 1, 0, label = 'R110')
# tran.add_fit(0, -2, 2, 0, label = 'R220')
# tran.add_fit(0, -3, 3, 0, label = 'R330')
# tran.add_fit(0, -4, 4, 0, label = 'R440')


tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

theta = np.pi/4
Ldip = rho*theta
beta1, beta2 = 0, 0

tran.add_rotation_beam(180)
tran.add_rotation(beta2, 'R21')
#tran.add_bend(0.62831, 60, B1, 'B1')
tran.add_bend(Ldip, theta/np.pi*180, 0, 'B2')
tran.add_rotation(beta1, 'R22')
tran.add_rotation_beam(-180)


tran.add_drift(D21, 'D21')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_quad(Lquad, Q21, Rquad, 'Q21')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D22, 'D22')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_quad(Lquad, Q22, Rquad, 'Q22')

tran.add_text()
tran.add_fit(2, 1, 1, xmax, label = 'F112')
tran.add_fit(2, 3, 3, ymax, label = 'F332')
tran.add_text()

tran.add_drift(D22, 'D22')

tran.add_quad(Lquad, Q23, Rquad, 'Q23')

tran.add_drift(D24, 'D21')


theta = np.pi/2
Ldip = rho*theta
beta1, beta2 = 0, 0

tran.add_rotation(beta2, 'R31')
#tran.add_bend(0.62831, 60, B1, 'B1')
tran.add_bend(Ldip, theta/np.pi*180, 0, 'B3')
tran.add_rotation(beta1, 'R32')


#tran.add_text()
#tran.add_fit(2, 1, 1, xmax, label = 'FIT1')
#tran.add_fit(2, 3, 3, ymax, label = 'FIT2')
#tran.add_text()

tran.add_fit(0, -1, 6, 0, label = 'R160')
tran.add_fit(0, -2, 6, 0, label = 'R260')

#tran.add_fit(0, 1, 1, xmatch, label = 'FIT3')
#tran.add_fit(0, 3, 3, ymatch, label = 'FIT4')
#tran.add_fit(2, 5, 5, zmatch, label = 'FIT5')
#tran.add_text()

tran.add_fit(0, 1, 1, xmatch, label = 'F110')
tran.add_fit(0, 3, 3, ymatch, label = 'F330')
tran.add_fit(0, 5, 5, zmatch, label = 'F550')
tran.add_text()

for _ in np.arange(5):
    tran.add_drift(0.1, 'D')
    tran.add_fit(0, -1, 6, 0, label = 'FI16')
    tran.add_fit(0, -2, 6, 0, label = 'R260')
    tran.add_text()
    
    tran.add_fit(0, 1, 1, xmatch, label = 'F110')
    tran.add_fit(0, 3, 3, ymatch, label = 'F330')
    tran.add_fit(0, 5, 5, zmatch, label = 'F550')
    tran.add_text()


# tran.add_quad(Lquad, Q11, Rquad, 'Q11')
# tran.add_drift(D12, 'D12')
# #tran.add_fit(0, -1, 6, 0, label = 'FI16')
# tran.add_fit(0, 1, 1, xmatch, label = 'FIT3')
# tran.add_fit(0, 3, 3, ymatch, label = 'FIT4')
# tran.add_fit(2, 5, 5, zmatch, label = 'FIT5')
# tran.add_text()

# tran.add_drift(D12, 'D12')
# tran.add_quad(Lquad, Q12, Rquad, 'Q12')
# #tran.add_fit(0, -1, 6, 0, label = 'FI16')
# tran.add_fit(0, 1, 1, xmatch, label = 'FIT3')
# tran.add_fit(0, 3, 3, ymatch, label = 'FIT4')
# tran.add_fit(2, 5, 5, zmatch, label = 'FIT5')
# tran.add_text()

# tran.add_drift(0.1, 'D13')
# #tran.add_fit(0, -1, 6, 0, label = 'FI16')
# tran.add_fit(0, 1, 1, xmatch, label = 'FIT3')
# tran.add_fit(0, 3, 3, ymatch, label = 'FIT4')
# tran.add_fit(2, 5, 5, zmatch, label = 'FIT5')
# tran.add_text()


tran.add_fit_emittance(2, xemi)
tran.add_sentinel()

tran.add_vary_code(5, '0F', 'Q11')
tran.add_vary_code(5, '0G', 'Q12')
tran.add_vary_code(5, '0F', 'Q13')
tran.add_vary_code(5, '0L', 'Q21')
tran.add_vary_code(5, '0M', 'Q22')
tran.add_vary_code(5, '0L', 'Q23')

tran.add_vary_code(3, '-B', 'D11')
tran.add_vary_code(3, 'B', 'D12')
#tran.add_vary_code(3, '-C', 'D21')
#tran.add_vary_code(3, 'C', 'D22')

#tran.add_vary_code(3, '-B', 'D1')
#tran.add_vary_code(3, 'B', 'D2')

# tran.add_vary_code(5,'0M','Q11')
# tran.add_vary_code(5,'0N','Q12')


tran.add_text()

tran.add_fitting(2, 'F112')
tran.add_fitting(2, 'F332')
tran.add_fitting(0, 'F110')
tran.add_fitting(0, 'F330')
tran.add_fitting(0, 'F550')

tran.add_fitting(0, 'R160')
tran.add_fitting(0, 'R260')

# tran.add_fitting(0, 'R110')
# tran.add_fitting(0, 'R220')
# tran.add_fitting(0, 'R330')
# tran.add_fitting(0, 'R440')

tran.add_fitting(2,'XEMI')
tran.add_text()

tran.add_sentinel('\nSENTINEL\nSENTINEL\n')

tran.output('C:\\TRANS\\for001.dat')
#os.system('cp for001.dat C:\\TRANS\\for001.dat')