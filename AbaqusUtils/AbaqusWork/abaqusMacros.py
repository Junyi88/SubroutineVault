# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

def Macro1():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    odb = session.odbs['Y:/Tabassam1/Creep2.odb']
    session.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
        variable=(('PEEQ', INTEGRATION_POINT), ('S', INTEGRATION_POINT), ), 
        elementPick=(('PART-1-1', 1, ('[#1 ]', )), ), )
    session.viewports['Viewport: 1'].view.setValues(session.views['User-4'])
    x0 = session.xyDataObjects['PEEQ PI: PART-1-1 E: 1 IP: 1']
    x1 = session.xyDataObjects['S:Max Principal (Abs) PI: PART-1-1 E: 1 IP: 1']
    x2 = session.xyDataObjects['S:Max Principal PI: PART-1-1 E: 1 IP: 1']
    x3 = session.xyDataObjects['S:Mid Principal PI: PART-1-1 E: 1 IP: 1']
    x4 = session.xyDataObjects['S:Min Principal PI: PART-1-1 E: 1 IP: 1']
    x5 = session.xyDataObjects['S:Mises PI: PART-1-1 E: 1 IP: 1']
    x6 = session.xyDataObjects['S:Pressure PI: PART-1-1 E: 1 IP: 1']
    x7 = session.xyDataObjects['S:S11 PI: PART-1-1 E: 1 IP: 1']
    x8 = session.xyDataObjects['S:S12 PI: PART-1-1 E: 1 IP: 1']
    x9 = session.xyDataObjects['S:S13 PI: PART-1-1 E: 1 IP: 1']
    x10 = session.xyDataObjects['S:S22 PI: PART-1-1 E: 1 IP: 1']
    x11 = session.xyDataObjects['S:S23 PI: PART-1-1 E: 1 IP: 1']
    x12 = session.xyDataObjects['S:S33 PI: PART-1-1 E: 1 IP: 1']
    x13 = session.xyDataObjects['S:Third Invariant PI: PART-1-1 E: 1 IP: 1']
    x14 = session.xyDataObjects['S:Tresca PI: PART-1-1 E: 1 IP: 1']
    session.writeXYReport(fileName='abaqus.rpt', xyData=(x0, x1, x2, x3, x4, x5, 
        x6, x7, x8, x9, x10, x11, x12, x13, x14))


def Macro2():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    odb = session.odbs['Y:/GitKrakenRes/AbaqusUtils/AbaqusWork/Exclude/MgO_10MicronC.odb']
    session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
        NODAL), ), nodeSets=('PART1', ))


