import wx,sys,os
import numpy as np
import raypath as rp
import laser as las

TITLE = "RayView"

class OptiTest(wx.Panel):
    def __init__(self,parent=None,size=(800,600)):
        wx.Panel.__init__(self,parent=parent,id=-1,size=size)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.opti = rp.Opti()
        self.optiDict = {'xPos':300,'nOpti':1,'Diameter':100,
                         'lFocal':200,'rFocal':200}
        self.sourceDict = {'xSource':10,'ySource':100,'theta':0}

        self.source = (20,30)
        self.abcs = []
        self.nodes = []
        self.dots = []
        self.rays = []
        self.InitPanel()

    def InitPanel(self):
        self.drawPanel = wx.Panel(self)

        #########初始化paraBook
        paraBook = wx.Notebook(self,size=(300,-1))
        optiPanel = wx.Panel(paraBook)
        sourcePanel = wx.Panel(paraBook)
        paraBook.AddPage(optiPanel,'opti')
        paraBook.AddPage(sourcePanel,'source')

        ###需要初始化edge
        self.setEdge()

        ####################optiBox###################
        self.paraSliders = {}
        optiBox = wx.BoxSizer(wx.VERTICAL)
        for key in self.optiDict:
            self.paraSliders[key]=wx.Slider(
                optiPanel,minValue=1,maxValue=1000,size=(200,-1))
            self.paraSliders[key].Bind(wx.EVT_SCROLL,
                lambda evt,mark=key: self.OnSliderScroll(evt,mark))
            optiBox.Add(self.paraSliders[key],proportion=1,
                flag=wx.LEFT|wx.CENTER)
            optiBox.Add(wx.StaticText(optiPanel,size=(120,30),label=key,
                style=wx.ALIGN_RIGHT),proportion=1,
                flag=wx.ALIGN_CENTER, border=10)

        self.testFlag = wx.TextCtrl(
            optiPanel,size=(250,400),value='hellos',style=wx.TE_MULTILINE)
        optiBox.Add(self.testFlag,proportion=1,
            flag=wx.ALIGN_CENTER|wx.ALL|wx.ALIGN_RIGHT,border=0)

        optiPanel.SetSizer(optiBox)

        ####################sourceBox###################

        sourceBox = wx.BoxSizer(wx.VERTICAL)
        for key in self.sourceDict:
            self.paraSliders[key]=wx.Slider(
                sourcePanel,minValue=1,maxValue=1000,size=(200,-1))
            self.paraSliders[key].Bind(wx.EVT_SCROLL,
                lambda evt,mark=key: self.OnSliderScroll(evt,mark))
            sourceBox.Add(self.paraSliders[key],proportion=1,
                flag=wx.LEFT|wx.CENTER)
            sourceBox.Add(wx.StaticText(sourcePanel,size=(120,30),label=key,
                style=wx.ALIGN_RIGHT),proportion=1,
                flag=wx.ALIGN_CENTER, border=10)

        sourcePanel.SetSizer(sourceBox)

        mainBox = wx.BoxSizer()
        mainBox.Add(self.drawPanel,proportion=1,flag=wx.ALL|wx.EXPAND,border=10)
        mainBox.Add(paraBook,proportion=0,flag=wx.ALL|wx.EXPAND,border=10)
        self.SetSizer(mainBox)



    def OnSliderScroll(self,evt,mark):
        paraArea = {'ySource':[-300,300],'xSource':[0,1000],
                    'xPos':[0,1200],'Diameter':[0,500],
                    'lFocal':[-1000,1000],'rFocal':[-1000,1000],
                    'theta':[0,np.pi*2],'nOpti':[0.1,10]}
        pValue = self.paraSliders[mark].GetValue()
        pMin,pMax=paraArea[mark]
        if mark in self.optiDict:
            self.optiDict[mark] = pMin+(pMax-pMin)/1000*pValue
        elif mark in self.sourceDict:
            self.sourceDict[mark]=pMin+(pMax-pMin)/1000*pValue
        pStr = ''
        for key in self.optiDict:
            pStr += key+':'+str(self.optiDict[key])+'\n'

        self.setEdge()
        self.getRay()
        pStr += '-----edge----------\n' + str(self.opti.edges) + '\n'

        pStr += '-----rays------\n'
        for ray in self.rays:
            pStr += str(ray)+'\n'

        pStr += '-----abcs------\n'
        for abc in self.abcs:
            pStr += str(abc)+'\n'

        self.testFlag.SetValue(pStr)
        self.DrawPath()

    def setEdge(self):
        lFocal = self.optiDict['lFocal']

        rFocal = self.optiDict['rFocal']

        self.opti.asLens([lFocal,rFocal],self.optiDict['Diameter'])
        self.opti.doMove((self.optiDict['xPos'],0))
        self.opti.setND(self.optiDict['nOpti'])
        self.source = (self.sourceDict['xSource'],
                       self.sourceDict['ySource'])

    def setSource(self):
        pass

    def getRay(self):
        self.rays,self.abcs,self.dots = [[],[],[]]

        sDot = self.source
        sRay = rp.getABC(self.sourceDict['theta'],sDot)
        inPoint,outPoint,flec,frac = self.opti.singleReflect(sRay,sDot,1)

        if inPoint == []: return [] 
        self.dots.append(inPoint)
        self.rays.append([sDot,inPoint])

        crossflec = self.crossRagion(flec,inPoint)
        if crossflec != []:
            self.dots.append(crossflec)
            self.rays.append([inPoint,crossflec])
            self.abcs.append(flec)

        if outPoint == []: return []
        self.dots.append(outPoint)
        self.rays.append([inPoint,outPoint])

        if frac == []: return []
        crossfrac = self.crossRagion(frac,outPoint)
        if crossflec != []:
            self.dots.append(crossfrac)
            self.rays.append([outPoint,crossfrac])
            self.abcs.append(frac)

    def crossRagion(self,ray,point):
        w,h = self.drawPanel.GetSize()
        edges = [[(0,0),(0,w)],[(0,h/2),(0,-h/2)],[(0,-h/2),(w,-h/2)],
                [(w,-h/2),(w,h/2)],[(w,h/2),(0,h/2)]]
        for dots in edges:
            cross=rp.getCross(ray,dots,point)
            if cross!=[]:
                return cross
        return []

    def OnPaint(self,evt):
        self.DrawPath() 

    def DrawPath(self):
        w,h = self.drawPanel.GetSize()

        dc = wx.ClientDC(self.drawPanel)
        dc.SetPen(wx.Pen('#666666'))
        dc.DrawRectangle(0,0,w,h)
        dc.SetDeviceOrigin(0,h/2)
        dc.SetAxisOrientation(True,True)

        dc.SetPen(wx.Pen('#0000FF'))
        dc.DrawLine(0,0,w,0)

        dc.SetPen(wx.Pen('#00FF00'))
        for edge in self.opti.edges:
            dots = edge['dots']
            if len(dots)==2:
                dc.DrawLine(dots[0][0],dots[0][1],
                            dots[1][0],dots[1][1])
            elif len(dots)==3:
                x3,y3,_=rp.arc2cir(dots)
                if dots[1][0]>dots[2][0]:
                    dc.DrawArc(dots[0][0],dots[0][1],
                               dots[1][0],dots[1][1],x3,y3)
                else:
                    dc.DrawArc(dots[1][0],dots[1][1],
                               dots[0][0],dots[0][1],x3,y3)


        dc.SetPen(wx.Pen('#FF0000'))
        dc.DrawCircle(self.source[0],self.source[1],10)
        for ray in self.rays:
            dc.DrawLine(ray[0][0],ray[0][1],
                        ray[1][0],ray[1][1])

        dc.SetPen(wx.Pen('#FF00FF'))
        for dot in self.dots:
            dc.DrawCircle(dot[0],dot[1],5)

#光路Frame
class RayFrame(wx.Frame):
    def __init__(self,Parent=None,title=TITLE):
        wx.Frame.__init__(self,parent=Parent,title=title,size=(800,600))
        self.mode = 'test'

        self.InitMode()
        self.Center()
        #self.Bind(wx.EVT_SIZE,self.OnSize)
        #self.Bind(wx.EVT_MOVE,self.OnSize)
        self.InitMenu()

    def InitMode(self):
        self.RayPanel = OptiTest(self)

    def InitMenu(self):
        self.mb = wx.MenuBar()

        ######fileMenu#######
        fileMenu = wx.Menu()
        fileMenu.Append(wx.ID_NEW,'New')
        fileMenu.Append(wx.ID_OPEN,'Open')
        fileMenu.Append(wx.ID_SAVE,'Save')
        fileMenu.AppendSeparator()
        fileMenu.Append(wx.ID_EXIT,'Quit')
        self.mb.Append(fileMenu,'file')

        self.Bind(wx.EVT_MENU,self.OnNew,id=wx.ID_NEW)
        self.Bind(wx.EVT_MENU,self.OnOpen,id=wx.ID_OPEN)
        self.Bind(wx.EVT_MENU,self.OnSave,id=wx.ID_SAVE)
        self.Bind(wx.EVT_MENU,self.OnQuit,id=wx.ID_EXIT)

        ######insertMenu#######
        insertMenu = wx.Menu()
        ID_RAYSOURCE = wx.NewId()
        ID_RAGION = wx.NewId()
        ID_PANE = wx.NewId()
        insertMenu.Append(ID_RAYSOURCE,'RaySource')
        insertMenu.Append(ID_RAGION,'Ragion')
        insertMenu.Append(ID_PANE,'Pane')
        self.mb.Append(insertMenu,'insert')

        self.Bind(wx.EVT_MENU,self.OnSetRaySource,id=ID_RAYSOURCE)
        self.Bind(wx.EVT_MENU,self.OnSetRagion,id=ID_RAGION)
        self.Bind(wx.EVT_MENU,self.OnSetPane,id=ID_PANE)


        ######modeMenu#######
        modeMenu = wx.Menu()
        modes = ['edge','test','path','mesure','wave','cavity']
        modeID = {}
        for mode in modes:
            modeID[mode] = wx.NewId()
            modeMenu.Append(modeID[mode],mode,kind=wx.ITEM_RADIO)
            self.Bind(wx.EVT_MENU,
            lambda evt,mark=mode:self.OnModeChange(evt,mark),id=modeID[mode])
        self.mb.Append(modeMenu,'mode')

        self.SetMenuBar(self.mb)

    def OnNew(self,evt):
        pass

    def OnOpen(self,evt):
        filesFilter = "opath(*.opath)|*.opath|""All files(*.*)|*.*"
        fileDialog = wx.FileDialog(self,message="select a file", wildcard=filesFilter,style=wx.FD_OPEN)
        dialogResult = fileDialog.ShowModal()
        fileName = fileDialog.GetFilename()
        dirName =  fileDialog.GetDirectory()
        if dialogResult != wx.ID_OPEN:
            filePath = dirName + '/' + fileName
            self.openRayPanel(filePath)

    def OnSave(self,evt):
        filesFilter = "opath(*.opath)|*.opath|""All files(*.*)|*.*"
        fileDialog = wx.FileDialog(self,message="select your saving path", wildcard=filesFilter,style=wx.FD_SAVE)
        dialogResult = fileDialog.ShowModal()
        fileName = fileDialog.GetFilename()
        dirName =  fileDialog.GetDirectory()
        if dialogResult != wx.ID_SAVE:
            filePath = dirName + '/' + fileName
            self.saveRayPanel(filePath)

    def openRayPanel(self,filePath):
        file = open(filePath,'r')
        sDict = eval(file.read())
        file.close()
        self.RayPanel.Destroy()
        self.mode = sDict['mode']
        self.RayPanel.optis = sDict['optis']

    def saveRayPanel(self,filePath):
        file = open(filePath,'w')
        sDict = {'mode':self.mode,'optis':self.RayPanel.optis}
        file.write(str(sDict))
        file.close()


        if self.mode=='test':
            pass

    def OnQuit(self,evt):
        pass

    def OnSetRaySource(self,evt):
        pass

    def OnSetRagion(self,evt):
        pass

    def OnSetPane(self,evt):
        pass

    def OnModeChange(self,evt,mark):
        self.RayPanel.Destroy()
        self.mode=mark
        self.InitMode()


    def OnSize(self,evt):
        self.RayPanel.absPos=self.Position
        self.Refresh()

class MyApp(wx.App):
    def OnInit(self):
        self.SetAppName(TITLE)
        self.Frame = RayFrame(None)
        self.Frame.Show()
        return True

if __name__ == '__main__':
    app = MyApp()
    app.MainLoop()