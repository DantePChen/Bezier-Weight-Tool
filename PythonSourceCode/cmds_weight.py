import maya.cmds as cmds
from maya.OpenMaya import *
from maya.OpenMayaAnim import *
import BezierWeightToolCompute as C_compute


def distance_between_joints(joint1,joint2):
    translation1 = cmds.xform(joint1, query = True, translation = True, worldSpace = True)
    translation2 = cmds.xform(joint2, query = True, translation = True, worldSpace = True)
    distance = ((translation2[0] - translation1[0])**2 + (translation2[1] - translation1[1])**2 + (translation2[2] - translation1[2])**2) ** 0.5
    return distance

def get_Geometry(geometry = None, shape_type = "mesh",intermediate=False):
    if geometry is None:
        selected = cmds.ls(selection=True)
        if len(selected)==0:
            return cmds.warning("please select a" + shape_type)
        geometry = selected[0]

    if cmds.nodeType(geometry)=="transform":
        shapes=cmds.listRelatives(geometry,shapes=True,path=True)
        if not shapes:
            return cmds.warning("there is no " + shape_type)
        for shape in shapes:
            is_intermediate = cmds.getAttr('{}.intermediateObject'.format(shape))
            if intermediate and is_intermediate and cmds.listConnections('{}.worldMesh'.format(shape), source=False):
                return shape
            elif not intermediate and not is_intermediate:
                return shape
    elif cmds.nodeType(geometry) in ["mesh","nurbsCurve","nurbsSurface"]:
        return geometry


def get_skinCluster(polygon=None):
    # first use function to get polygon
    if polygon is None:
        polygon = get_Geometry(shape_type="mesh")

    # then search from history to find skinCluster
    if polygon is None:
        return

    history_nodes = cmds.listHistory(polygon)
    for history in history_nodes:

        if cmds.nodeType(history) == "skinCluster":
            return history

    cmds.warning("Can't find skinCluster")

def setToVertexMode():
    cmds.setToolTo('selectSuperContext')
    obj = cmds.ls(selection=True)
    cmds.selectType(objectComponent=True, allComponents=False)
    cmds.selectType(objectComponent=True, vertex=True)
    cmds.selectType(vertex=True)
    cmds.hilite(obj)

def get_single_joint_info():
    polygon = get_Geometry()
    sk = get_skinCluster(polygon)
    """
    liw is Limit Information Weight 
    if liw is 1, means its weight can be painted,
    so when it's 0, means the joint is locked
    """
    unlock_joints = []
    influences = cmds.skinCluster(sk,query = True, influence = True)
    for joint in influences:
        if cmds.nodeType(joint) == "joint" and not cmds.getAttr("{}.liw".format(joint)):
            unlock_joints.append(joint)
    if len(unlock_joints) != 1:
        return cmds.warning("You can have only one unlock joint!")
    unlock_joint = unlock_joints[0]

    paint_joints = []
    input_nodes = cmds.listConnections("{}.{}".format(sk, "paintTrans"), source=True, destination=False)
    for joint in input_nodes:
        if cmds.nodeType(joint) == "joint":
            paint_joints.append(joint)
    if len(paint_joints) != 1:
        return cmds.warning("You can have only one paint joint!")
    paint_joint = paint_joints[0]
    # get the joint world matrix
    paint_joint_M = cmds.xform(paint_joint, query=True, matrix=True, worldSpace=True)
    # get all the points from shape
    pointList = cmds.xform("{}.vtx[*]".format(polygon), query=True, translation=True, worldSpace=True)
    return [polygon,sk,paint_joint,paint_joint_M,pointList]

def get_paint_joint():
    polygon = get_Geometry()
    sk = get_skinCluster(polygon)
    paint_joints = []
    input_nodes = cmds.listConnections("{}.{}".format(sk, "paintTrans"), source=True, destination=False)
    for joint in input_nodes:
        if cmds.nodeType(joint) == "joint":
            paint_joints.append(joint)
    paint_joint = paint_joints[0]
    return paint_joint

def get_joint_index(skin_cluster,joint):
    influence_objects = cmds.skinCluster(skin_cluster, query=True, influence=True)

    try:
        joint_index = influence_objects.index(joint)
        return joint_index
    except ValueError:
        cmds.warning("Joint '{}' not found in the influence objects of the skin cluster.".format(joint))
        return None

def paint_weight(polygon,skin_cluster,paint_joint,weights):
    selected = MSelectionList()
    selected.add(skin_cluster)
    selected.add(polygon + ".vtx[*]")
    depend_node =MObject()
    selected.getDependNode(0,depend_node)
    fn_skin = MFnSkinCluster(depend_node)
    dag_path = MDagPath()
    components = MObject()
    api_weights= MDoubleArray()
    api_weights.setLength(len(weights))
    for i,w in enumerate(weights):
        api_weights[i] = weights[i]
    index = get_joint_index(skin_cluster,paint_joint)
    joint_index = MIntArray()
    joint_index.setLength(1)
    joint_index[0]=index
    selected.getDagPath(1,dag_path,components)
    fn_skin.setWeights(dag_path,components,joint_index,api_weights,True)
    refresh_weight(paint_joint)

def refresh_weight(paint_joint):
    selected_object = cmds.ls(selection=True)
    cmds.select(paint_joint)
    cmds.select(selected_object[0])

def paint_SingleWeights(px,py,radius,axis,infoList):
    polygon = infoList[0]
    skin_cluster = infoList[1]
    paint_joint = infoList[2]
    paint_joint_M = infoList[3]
    pointList =infoList[4]

    localPx = C_compute.WorldTransLocal(paint_joint_M,pointList,axis)

    if axis==3:
        weights = C_compute.get_Weights(localPx, px, py, radius, 2)
    else:
        weights = C_compute.get_Weights(localPx,px,py,radius,1)

    paint_weight(polygon,skin_cluster,paint_joint,weights)


def get_soft_localPx(radius=None):
    polygon = get_Geometry()
    skin_cluster = get_skinCluster(polygon)
    polygon = polygon.split(".")[0]
    cmds.softSelect(sse=1)
    cmds.softSelect(ssc="0,1,1,1,0,1")
    if radius is None:
        radius = cmds.softSelect(q=1, ssd=1)

    cmds.softSelect(ssd=radius)
    old_points = cmds.xform("{}.vtx[*]".format(polygon), query=True, translation=True, worldSpace=True)
    cmds.move(1,y=True,relative=True)

    new_points = cmds.xform("{}.vtx[*]".format(polygon), query=True, translation=True, worldSpace=True)
    cmds.move(-1,y=True,relative=True)

    localPx = C_compute.get_Soft_LocalPx(new_points,old_points)

    return [localPx,skin_cluster,polygon]

def paint_soft(px,py,radius):
    if radius is None:
        templist = get_soft_localPx()
        radius = 1
    else:
        templist = get_soft_localPx(radius)
    localPx = templist[0]
    skin_cluster = templist[1]
    polygon = templist[2]
    weights = C_compute.get_Weights(localPx,px,py,radius,0)
    paint_joint = get_paint_joint()
    cmds.select(polygon)
    cmds.ArtPaintSkinWeightsTool()
    paint_weight(polygon,skin_cluster,paint_joint,weights)

def get_chain_info():
    polygon = get_Geometry(shape_type="mesh")
    skin_cluster = get_skinCluster(polygon)

    # get the order of joints chain
    unlock_joints = []
    influences = cmds.skinCluster(skin_cluster, query=True, influence=True)
    for joint in influences:
        if cmds.nodeType(joint) == "joint" and not cmds.getAttr("{}.liw".format(joint)):
            unlock_joints.append(joint)
    num_unlock_joints = len(unlock_joints)
    # create a dict, key is their distance, the value is the pair of joints
    if num_unlock_joints >1:
        distance_joints = {}
        for i, joint1 in enumerate(unlock_joints):
            for joint2 in unlock_joints[i + 1:]:
                distance = distance_between_joints(joint1,joint2)
                distance_joints[distance] = [joint1, joint2]
        # find the most close pair of joints
        sorted_joints = distance_joints[min(distance_joints.keys())]
        unlock_joints.remove(sorted_joints[0])
        unlock_joints.remove(sorted_joints[1])
        # construct the right order of joints chain
        # first create a dict to get the distance between the most close joints and other joints
        for i in range(len(unlock_joints)):
            distance ={distance_between_joints(joint1, joint2): [joint1, joint2]
                      for joint1 in unlock_joints for joint2 in [sorted_joints[0], sorted_joints[-1]]}
            # find the most close joints pair
            sorted_joint, start_or_end_joint = distance[min(distance.keys())]
            # if star or end joint is the first joint in sorted joints, put it in the start of sorted joints list
            if start_or_end_joint == sorted_joints[0]:
                sorted_joints.insert(0,sorted_joint)
            # or put it in the end of the sorted joints list
            else:
                sorted_joints.append(sorted_joint)
            # then remove the sorted joint from original joint list
            unlock_joints.remove(sorted_joint)
        # create an array get the joints index
        joints_index = MIntArray()
        joints_num = len(sorted_joints)
        joints_index.setLength(joints_num)
        for i, joint in enumerate(sorted_joints):
            joints_index[i] = influences.index(joint)
    else:
        sorted_joints = unlock_joints
        joints_index = MIntArray()
        joints_index[0] = influences.index(sorted_joints[0])
        joints_num = 1

    # get the original weights before start painting
    pointList = cmds.xform("{}.vtx[*]".format(polygon), query=True, translation=True, worldSpace=True)
    selected = MSelectionList()
    selected.add(skin_cluster)
    selected.add(polygon + ".vtx[*]")
    depend_node = MObject()
    selected.getDependNode(0, depend_node)
    fn_skin = MFnSkinCluster(depend_node)
    path = MDagPath()
    components = MObject()
    weights = MDoubleArray()
    selected.getDagPath(1, path, components)
    fn_skin.getWeights(path,components,joints_index,weights)
    originalWeights_list = [weights[i] for i in range(weights.length())]
    #print(originalWeights_list)

    # get the joints matrix
    joints_Matrix=[]
    for joint in sorted_joints:
        paint_joint_M = cmds.xform(joint, query=True, matrix=True, worldSpace=True)
        joints_Matrix+=paint_joint_M

    return [pointList,originalWeights_list,joints_Matrix,sorted_joints,joints_index,components,path,fn_skin]

def paint_Chain_weights(px,py,radius,axis,ChainInfoList):
    pointlist = ChainInfoList[0]
    originalWeight_List = ChainInfoList[1]
    joints_Matrix =ChainInfoList[2]
    unlock_joints = ChainInfoList[3]
    joints_index = ChainInfoList[4]
    components = ChainInfoList[5]

    path = ChainInfoList[6]

    fn_skin = ChainInfoList[7]
    c_reault_weight = C_compute.get_Chain_Weights(joints_Matrix,pointlist,originalWeight_List,px,py,radius,1,axis)
    api_weights = MDoubleArray()
    api_weights.setLength(len(c_reault_weight))
    for i,w in enumerate(c_reault_weight):
        api_weights[i] = c_reault_weight[i]
    fn_skin.setWeights(path,components,joints_index,api_weights,1)
    refresh_weight(unlock_joints[0])