#coding=utf-8
import os
import maya.cmds as cmds
import maya.mel as mel

def shelfButtonInstall():
    # get the absolute path of current file
    Path = os.path.dirname(__file__)
    # modify the path from system
    Path = os.path.abspath(Path).replace("\\","/")
    Label = "BWT"
    Script = '''
import sys
if \"%s\" not in sys.path:
    sys.path.append(\"%s\")
    
import %s as UI

reload(UI)

UI.show() 
''' % (Path,Path,"BezierUI")

    mel.eval('global string $gShelfTopLevel')
    gShelfTopLevel = mel.eval('$tmp = $gShelfTopLevel')

    currentShelf = cmds.tabLayout(gShelfTopLevel, query=True, selectTab=True)
    cmds.setParent(currentShelf)

    iconExt = "png"
    icon = "pythonFamily." + iconExt

    cmds.shelfButton(
        command=Script,
        annotation=Label,
        label=Label,
        imageOverlayLabel=Label,
        image=icon,
        image1=icon,
        sourceType="python"
    )
    print("Button created")

def onMayaDroppedPythonFile(param):
    shelfButtonInstall()