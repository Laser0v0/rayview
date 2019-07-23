import numpy as np 
import cmath
import sys

#移动dMove,dMove建议为元组
doMove = lambda dots,dMove : [
    (p[0]+dMove[0],p[1]+dMove[1]) for p in dots]

#绕中心旋转,flag=True时，第一个点是中心，否则无中心
def rotByCenter(dots,theta,center=[]):
    if center == []:
        center = np.array(dots[0])
    else:
        center = np.array(center)
    rotMat = np.array([[np.cos(theta),np.sin(theta)],
                       [-np.sin(theta),np.cos(theta)]])
    dots = [tuple(p+np.matmul((p-center),rotMat)) for p in dots]
    return dots

def doResize(dots,aSize,bSize,flag='single'):
    if flag=='scale':
        scale = np.array([aSize,bSize])
    else:
        scale = np.array(bSize).astype(float)/aSize
    if flag=='single':
        scale = min(scale)*np.array([1,1])
    dots = [tuple(p*scale) for p in dots]
    return dots

#[thetax=costheta,thetay=sintheta]
#取值与象限关系[[+,+],[-,+],[-,-],[+,-]]
def getTheta(abc=[1,1,1]):
    r = np.sqrt(abc[0]**2+abc[1]**2)
    return [-abc[1]/r,abc[0]/r]

#通过[thesaX,thetaY]
def getAngle(theta):
    theta = np.array(theta)/np.sqrt(theta[0]**2+theta[1]**2)
    if theta[1]>=0:
        return np.arccos(theta[0])
    else:
        return (np.pi-np.arccos(theta[0]))

#通过点对求角度
def getDotsAngle(dots=[(0,0),(1,1)]):
    return getAngle(np.array(dots[1])-dots[0])

#当输入一个值时，为点对#第一个点为始点，第二个点为终点
#当输入两个值时，为角度和始点theta=0,sPoint=(1,1)
##ab取值与象限关系[[+,-],[+,+],[-,+],[-,-]]
def getABC(*par):
    if len(par)==1:
        dots = par[0]
        abc = [dots[1][1]-dots[0][1],
                dots[0][0]-dots[1][0], 
                -np.linalg.det(dots)]
        return np.array(abc)/(np.sqrt(abc[0]**2+abc[1]**2))
    elif len(par)==2:
        theta,sPoint = par
        a,b = [np.sin(theta),-np.cos(theta)]
        c = -(a*sPoint[0]+b*sPoint[1])
        return [a,b,c]

#将点按照射线的方向排序
def sortByRay(abc,points):
    ab = np.array([abc[1],-abc[0]])
    pDict = {ab.dot(point):point for point in points}
    keyList = list(pDict.keys())
    keyList.sort(reverse=True)
    return [pDict[key] for key in keyList]

#直线和点集是否存在交点#isCross(abc,arc)
def isCross(abc,dots):
    abc = np.array(abc)
    flags = [abc.dot(list(p)+[1]) for p in dots]
    flag1,flag2,flag3 = [0,0,0]
    for flag in flags:
        if flag > 0:
            flag1 += 1
        elif flag <0:
            flag2 += 1
        else:
            flag3 += 1
    return flag1*flag2+flag3

#返回为0时表示在直线上
isOnLine = lambda abc,x,y : np.array(abc).dot([x,y,0])
#返回为1时，表示在射线上
isOnRay = lambda abc,points : sortByRay(abc,points)==points

#得到点关于直线的对称点,k为比例系数
def getSymDot(point,line,k=1):
    return tuple((np.array(getPedal(point,line))*(1+k)-point)/k)

#得到直线的垂足
def getPedal(point,line):
    a,b,c=line
    x0,y0 = point
    y1 = (a**2*y0-a*b*x0-b*c)/(a**2+b**2)
    x1 = (b**2*x0-a*b*y0-a*c)/(a**2+b**2)
    return (x1,y1)

#是否为劣弧#目前所有的弧都假定为劣弧，否则需要此判定函数
def isMinorArc(arc):
    arc = np.array(arc)
    flag = np.sum(arc[0]-arc[2])*(arc[1]-arc[2])
    if flag==0:
        return 0
    else:
        return -1 if flag <  0 else 1

#获取圆上的切线#arc2cir(arc)
def getTangent(corss=[0,1],circle=[0,0,1]):
    a = corss[0]-circle[0]
    b = corss[1]-circle[1]
    c = -a*corss[0]-b*corss[1]
    return [a,b,c]

#对于与cir的交点，是否在劣弧上
def isOnArc(point,arc):
    arc = np.array(arc)
    dAB = 0.5*np.linalg.norm(arc[0]-arc[1])
    dCrossA = np.linalg.norm(0.5*(arc[0]+arc[1])-point)
    return dAB > dCrossA

#弧转圆
def arc2cir(arc):
    arc = np.array(arc)
    dCD = np.linalg.norm(1/2*(arc[0]+arc[1])-arc[2])
    dBC2 = np.sum(np.square(arc[1]-arc[2]))
    radius = 0.5*dBC2/dCD
    theta = (arc[2]-1/2*(arc[0]+arc[1]))/dCD
    zero = arc[2]-radius*theta
    return list(zero)+[radius]

#直线和圆的交点
#输入abc,cir为list,point为tuple
#输出：无交点时输出p[],有交点时输出[()]
def getCrossCircle(abc=[1,1,1],cir=[0,0,2],point=()):
    c=np.array(cir[0:2]+[1]).dot(abc)
    a2b2 = abc[0]**2+abc[1]**2
    delt = a2b2*cir[2]**2-c**2
    if delt<0: return []
    else: delt=np.sqrt(delt)
    
    plusMinus = lambda a,b : np.array([a-b,a+b])
    yCross = plusMinus(-abc[1]*c,abc[0]*delt)/a2b2*[1,1]+cir[1]
    xCross = plusMinus(-abc[0]*c,-abc[1]*delt)/a2b2*[1,1]+cir[0]
    if point==[]:
        return [(xCross[i],yCross[i]) for i in [0,1]]
    yFlag = (yCross-point[1])*abc[0] >= 0
    xFlag = (point[0]-xCross)*abc[1] >= 0
    zFlag = np.abs(xFlag)+np.abs(yFlag) > 0
    flag = yFlag*xFlag*zFlag

    return [(xCross[i],yCross[i]) 
            for i in range(len(yFlag)) if flag[i]]

#the first cross between an arc and a ray from a point
def getCrossArc(abc=[1,1,1],arc=[[0,1],[0,-1],[1,0]],point=[]):
    if  point == []:
        return []
    crossDict = {np.sqrt((p[0]-point[0])**2+(p[1]-point[1])**2):p
                 for p in getCrossCircle(abc,arc2cir(arc),point) 
                 if (isOnArc(p,arc) and (p!=point))}
    if crossDict == {}:
        return []
    return  crossDict[min(crossDict.keys())]

# the cross between a ray and segment
# #不考虑端点#dots也可以为abc
def getCrossDots(abc=[1,-1,1],dots=[[0,2],[2,0]],point=[]):
    abc = np.array(abc)
    if len(dots)==2:
        flags = [abc.dot(list(p)+[1]) for p in dots]
        if flags[0]*flags[1]>=0: 
            return []       #此时无交点或者交于端点
        abc2 = getABC(dots)
    else: 
        abc2 = dots
    pVar= np.linalg.solve(
          np.array([abc[0:2],abc2[0:2]]),
          np.array([[-abc[2]],[-abc2[2]]]))     ##[[x],[y]]
    cross = (pVar[0][0],pVar[1][0])
    if point == []: 
        return cross                #此时为无端直线
    if cross == point:
        return []
    abcTest = getABC([point,cross])
    if (abcTest[0]*abc[0]>=0) and (abcTest[1]*abc[1]>=0):
        return cross
    return []

#abc and points
def getCross(abc=[1,-1,0],dots=[[0,-1],[0,1],[0.5,0]],point=[]):
    if len(dots)==3:
        return getCrossArc(abc,dots,point)
    if len(dots)==2:
        return getCrossDots(abc,dots,point)

#
def getCrossEdges(abc=[1,1,1],edges=[[(0,-1),(0,1)],[(0,1),(0,-1),(1/2,0)]],point=[-1,0]):
    crosses = []
    crossEdges = []
    for edge in edges:
        cross = getCross(abc,edge,point)
        if cross != []:
            crosses.append(cross)
            crossEdges.append(edge)
    if crosses == []:
        return [],[]
    crossDict = {(abs(crosses[i][0]-point[0])+abs(crosses[i][0]-point[0])):i 
                  for i in range(len(crosses))}
    crossIndex = crossDict[min(crossDict.keys())]
    return crosses[crossIndex],crossEdges[crossIndex]

#光在直线表面的反折射
def cataDioLine(abc=[1,-1,1],line=[2,-1,1],
                sPoint=[],cross=[],n1=1,n2=1.5):
    normal = [-line[1],line[0],line[1]*cross[0]-line[0]*cross[1]]#法线
    flecDot = getSymDot(sPoint,normal)
    flec=getABC([cross,flecDot])

    dPara = np.sqrt(line[0]**2+line[1]**2)
    dNormal = np.abs(np.array(normal).dot(list(sPoint)+[1]))/dPara#到法线距离
    dPane = np.abs(np.array(line).dot(list(sPoint)+[1]))/dPara#到反射面距离

    if dNormal == 0:
        return flec,abc
    delt = (n2/n1)**2*(1+(dPane/dNormal)**2)-1#判定全反射
    if delt>0:
        k =dPane/dNormal/np.sqrt(delt)
        fracDot = getSymDot(sPoint,normal,k)
        fracDot = getSymDot(fracDot,line)
        frac = getABC([cross,fracDot])
        return flec,frac
    return flec,[0,0,0]

#光在直线或弧线表面的反折射
def cataDio(abc=[1,-1,1],dots=[(0,2),(2,2)],
            sPoint=[-2,-1],n1=1,n2=1.5):
    cross = getCross(abc,dots,sPoint)
    if cross == []:
        return [],[],[]
    if len(dots)==3:
        line = getTangent(cross,arc2cir(dots))
    elif len(dots)==2:
        line = getABC(dots)
    flec,frac = cataDioLine(abc,line,sPoint,cross,n1,n2)
    return cross,flec,frac

#输入为点对的列表[points1,points2,..],点对用元组表示
def getAllPoints(edges=[((0,-1),(0,1)),((0,1),(0,-1),(1/2,0))]):
    outDots = []
    for dots in edges:
        outDots += list(dots)
    return list(set(outDots))

#生成透镜
def setLens(f=[1,1],D=1,thick=0):
    r = 2*np.array([-1,1])*f
    y = D/2
    x = r*np.sqrt(1-(y/r)**2)
    edge1,edge2 = [[(0,y),(0,-y),(r[0]-x[0],0)],
                   [(0,y),(0,-y),(r[1]-x[1],0)]]
    if (abs(r[0])<y):
        edge1 = [(0,y),(0,-y)]
    if (abs(r[1])<y):
        edge2 = [(0,y),(0,-y)]
    if thick==0:
        return [edge1,edge2]
    else:
        offset = thick/2
        edge1 = doMove(edge1,(-offset,0))
        edge2 = doMove(edge2,(offset,0))
        edge3 = [(-offset,y),(offset,y)]
        edge4 = [(-offset,-y),(offset,-y)]
        return [edge1,edge2,edge3,edge4]

#设置正多边形
def setRegular(num=5,d=10):
    r = d/np.sin(np.pi/num)
    dots = [(r*np.cos(i*2*np.pi/num),r*np.sin(i*2*np.pi/num))
            for i in range(num)]
    return [[dots[i%num],dots[(i+1)%num]] for i in range(num)]


#获取反射率和透过率,theta为弧度
def getFresnel(xTheta,yTheta, n1, n2):
    mid = np.sqrt(1-(n1/n2*yTheta)**2)#中间变量
    rp = (n2*xTheta-n1*mid)/(n2*xTheta+n1*mid)#p分量振幅反射率
    rs = (n1*xTheta-n2*mid)/(n1*xTheta+n2*mid)
    return rp**2,rs**2,1-rp**2,1-rs**2#分别是Rp, Rs, Tp, Ts

class Opti():
    def __init__(self,edges=[[(0,-1),(0,1)],[(0,1),(0,-1),(1/2,0)]],n=1.5):
        self.paras={'n':n,'centerless':True,'index':0}
        self.edges = [{'index':i,'dots':edges[i]} 
                      for i in range(len(edges))]
        self.rays = []
        self.setCenter()

    #edge格式为(dot1,dot2,...)
    def insertEdge(self,edge,albedo=0):
        self.edges.append(
            {'index':len(self.edges),'dots':edge})
        self.setCenter()

    #可接受编号和点集
    def delEdge(self,edge):
        try:
            if isinstance(edge,list):
                for edg in self.edges:
                    if edg['dots']==edge:
                        edge = edg['index']
            del self.edges[edge]
            self.update()
        except:
            print("no this edge")

    def doMove(self,pos):
        for i in range(len(self.edges)):
            self.edges[i]['dots'] = \
                doMove(self.edges[i]['dots'],pos)
        if not self.paras['centerless']:
            self.paras['center']=(self.paras['center'][0]+pos[0],
                                  self.paras['center'][1]+pos[1])
        else:
            self.setCenter()

    def doRot(self,theta):
        if self.paras['centerless']:
            self.setCenter()
        for i in range(self.edges):
            self.edges[i]['dots'] = \
                rotByCenter(self.edges[i]['dots'],theta,
                            self.paras['center'])

    def asLens(self,f=[1,1],D=1,thick=0):
        edges = setLens(f,D,thick)
        self.edges = [{'index':i,'dots':edges[i]} 
                      for i in range(len(edges))]
        self.setCenter((0,0))

    def setCenter(self,point=()):
        if point != []:
            self.paras['centerless'] = False
            self.paras['center']=point
        else:

            dots = getAllPoints([self.edges[i]['dots'] 
                                for i in range(len(self.edges))])
            center = np.mean(np.array(dots),0)
            self.paras['center'] = (center[0],center[1])
            self.paras['centerless']= True

    def setND(self,n):
        self.paras['n'] = n

    def setFilm(self,index=0,dWave=0,value=0):
        self.edges[index][dWave]=value

    #filmList格式：[[dWave,Value],...]
    def setFilms(self,filmList):
        for i in range(filmList):
            self.edges[i][filmList[0]]=filmList[1]

    def isCrossed(self,abc,point):
        pass

    def getCross(self,abc,point):
        dCross = 1e15
        index = False
        for edge in self.edges:
            cross = getCross(abc,edge['dots'],point)
            if cross != []:
                midDis = (point[0]-cross[0])**2+(point[1]-cross[1])**2
                if dCross > midDis:
                    dCross = midDis
                    index = edge['index']
        if not index: return [],[]
        return cross,self.edges[index]

    def singleReflect(self,abc,point,n):
        edges = [edge['dots'] for edge in self.edges]
        cross,edge0 = getCrossEdges(abc,edges,point)
        if cross==[]:
            return [],[],[],[]
        self.inPoint,self.flec,frac = cataDio(
            abc,edge0,point,n,self.paras['n'])
        cross,edge1 = getCrossEdges(abc,edges,self.inPoint)
        self.outPoint,_,self.frac = cataDio(
            frac,edge1,self.inPoint,self.paras['n'],n)
        return self.inPoint,self.outPoint,self.flec,self.frac

    def multiReflect(self,abc,point,n):
        #edge = self.getCross(abc,point)
        pass

def raySource(pos=(1,1),mode='single',para=[0],abc=(1,1,1)):
    if mode=='single':
        abc = getABC(para[0],pos)
        return [{'points':[pos],'theta': para[0],'abc':abc}]
    if mode=='pane':
        #此时的输入参数格式为[角度,孔径,直线数]
        theta = para[0]
        ab = [np.sin(theta),-np.cos(theta)]
        dUnit = para[1]/para[2]*ab
        pos = pos-para[1]/2*ab
        return [{'points':[pos+dUnit*i],'theta': theta,
            'abc':ab+[-dUnit*i]} for i in range(para[2])]
    if mode=='shperical':
        #此时输入参数格式为[theta1,theta2,n]    
        aUnit = (para[1]-para[0])/(para[2]-1)
        return [{'points':[pos],'theta': para[0]+aUnit*i,
                 'abc':getABC(para[0]+aUnit*i,pos)} 
                 for i in range(para[2])]


class RayPath():
    def __init__(self):
        self.rays = []      #当前光线
        self.nodes = []     #光线与器件的交点
        self.edges = []     #光学曲面
        self.optis = []     #光学仪器，提供光学曲面和点的索引
        self.paras = {}     #系统参数

    #插入新曲面，edge可以为字典和点集
    def insertEdge(self,edge,pos=(0,0),theta=0):
        newEdge = {'index':len(self.edges),'opti':len(self.optis)}
        if isinstance(edge,dict):
            newEdge.update(edge)
        else:
            newEdge['dots']=edge

        if newEdge['opti']==len(self.optis):
            self.optis.append({'index':len(self.optis),'edges':[len(self.edges)]})
        else:
            self.optis[newEdge['opti']]['edges'].append(len(self.edges))

        self.edges.append(edge)

    def insertOpti(self,opti,pos,theta):
        opti.doRot(theta)
        opti.doMove(pos)
        opti.paras['index']=self.paras['nOpti']
        self.updatePara()

    def delOpti(self,opIndex):
        del self.optis[opIndex]
        for index in range(opIndex,self.paras['nOpti']-1):
            self.optis[index].paras['index'] = index

    def updatePara(self):
        self.paras['nOpti'] = len(self.optis)
        self.paras['nRay'] = len(self.rays)
        self.paras['nNode'] = len(self.nodes)

    #paraDict 可输入power,polar,dWave等
    def insertRay(self,ray,point,paraDict={}):
        node = {'index':self.paras['nNode'],'childs':[],'parent':False,
                'edges':False,'rays':False,'pos':point}
        ray = {'index':self.paras['nRays'],'sPoint':node['index'],'abc':ray}
        node['rays']=[ray['index']]
        if isinstance(paraDict,dict):
            self.nodes.append(dict(node,**paraDict))
            self.rays.append(dict(ray,**paraDict))
        else:
            print('type error')
            self.nodes.append(node)
            self.rays.append(ray)

    def insertSource(self,ray):
        self.lives = ray
        pass

    def delSource(self):
        pass

    def oneSimpleStep(self):
        newRays = []
        for ray in self.rays:
            sNode = self.nodes[ray['node']]
            crossDict = {}
            for index in range(self.paras['nOpti']):
                if self.optis[index].isCrossed(ray,sNode['pos']):
                    cross = self.optis[index].getCross(ray,sNode['pos'])
                    crossDict[cross] = index
            opti = self.optis[sortByRay(ray,list(crossDict.keys()))[0]] #此为与该射线相交的第一个opti
            inPoint,outPoint,flec,frac = opti.simpleReflect(
                ray,sNode['pos'])
            if flec != False:
                self.updatePara()
                node = {'index':self.paras['nNode'],
                        'childs':[],'parent':sNode['index'],
                        'edges':opti['index'],'rays':[],
                        'pos':inPoint}
                self.nodes[sNode['index']]['childs'].append(node['index'])
                newRay = {'index':self.paras['nRays'],
                          'sPoint':node['index'],'abc':flec}
                node['rays'].append(ray['index'])
                self.nodes.append(node)
                newRays.append(newRay)
            if frac != False:
                self.updatePara()
                node = {'index':self.paras['nNode'],
                        'childs':[],'parent':sNode['index'],
                        'edges':opti['index'],'rays':[],
                        'pos':outPoint}
                self.nodes[sNode['index']]['childs'].append(node['index'])
                newRay = {'index':self.paras['nRays'],
                          'sPoint':node['index'],'abc':flec}
                node['rays'].append(ray['index'])
                self.nodes.append(node)
                newRays.append(newRay)
            self.nodes[ray['node']].remove([ray['index']])
        self.rays=newRay
        self.updatePara()

    def drawLines(self):
        lines = self.getLines(self.nodes[0])
        return [[self.nodes[ind[0]]['pos'],
                 self.nodes[ind[1]]['pos']] for ind in lines]

    def getLines(self,node):
        lines = []
        for index in node['childs']:
            lines += ([node['index'],index]) + self.getLines(self.nodes[index])
        return lines
