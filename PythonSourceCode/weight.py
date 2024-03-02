import pymel.core as pm
from maya.OpenMaya import *
from maya.OpenMayaAnim import *
import BezierWeightToolCompute as C_compute
import maya.cmds as cmds


def distance_between_joints(joint1,joint2):
    return (joint1.getTranslation(space="world")-joint2.getTranslation(space="world")).length()

def get_Geometry(geometry = None,shape_type = "mesh"):

    if geometry is None:
        selected = pm.selected(o=1)
        if len(selected) == 0:
            return pm.warning("please select a" + shape_type)
        geometry = selected[0]

    if geometry.type() == shape_type:
        return geometry.getParent()
    if geometry.type() != "transform":
        return pm.warning("please select a" + shape_type)
    shape = geometry.getShape()
    if not shape:
        return pm.warning("please select a" + shape_type)
    if shape.type() != shape_type:
        return pm.warning("please select a" + shape_type)
    return geometry

def get_skinCluster(polygon = None):
    # first use function to get polygon
    if polygon is None:
        polygon = get_Geometry(shape_type="mesh")
    # then search from history to find skinCluster
    if polygon is None:
        return
    for history in polygon.history(type = "skinCluster"):
        return history
    pm.warning("can't find skinCluster")


def refresh_weight():
    polygon = get_Geometry()
    sk = get_skinCluster(polygon)
    paint_joints = [joint for joint in sk.paintTrans.inputs() if joint.type() == "joint"]
    paint_joint = paint_joints[0]
    print(paint_joint)
    selcted_object = cmds.ls(selection=True)
    pm.select(paint_joint)
    cmds.select(selcted_object[0])

def update_paint():
    polygon = get_Geometry()
    sk = get_skinCluster(polygon)
    paint_joints = [joint for joint in sk.paintTrans.inputs() if joint.type() == "joint"]
    if len(paint_joints) != 1:
        return pm.mel.warning("\nyou need a paint joint")
    paint_joint = paint_joints[0]
    index = sk.influenceObjects().index(paint_joint)
    sk.ptw.set(sk.getWeights(polygon, index))

def get_Single_joint_info():
    polygon = get_Geometry()
    skin_cluster = get_skinCluster(polygon)
    """
    liw is Limit Information Weight 
    if liw is 1, means its weight can be painted,
    so when it's 0, means the joint is locked
    """
    unlock_joints = []
    for joint in skin_cluster.getInfluence():
        if joint.type()=="joint" and not joint.liw.get():
            unlock_joints.append(joint)
    if len(unlock_joints) !=1:
        return pm.mel.warning("You can have only one unlock joint!")
    unlock_joint = unlock_joints[0]

    paint_joints =[]
    for joint in skin_cluster.paintTrans.inputs():
        if joint.type()=="joint":
            paint_joints.append(joint)
    if len(paint_joints) !=1:
        return pm.mel.warning("You can have only one paint joint!")
    paint_joint = paint_joints[0]

    return paint_joint


def paint_SingleWeights(px,py,radius,axis):
    polygon = get_Geometry()
    skin_cluster = get_skinCluster(polygon)
    # get the paint joint
    paint_joint = get_Single_joint_info()

    # get the joint world matrix
    paint_joint_M = paint_joint.getMatrix(worldSpace=True)

    # convert it to a flat list for computing in C++
    flat_paint_joint_M = [float(value) for row in paint_joint_M for value in row]
    # get all the points from shape
    pointList = polygon.getShape().getPoints(space = "world")
    # convert it to a flat list for computing in C++
    flatFloatPointList = [float(coord) for point in pointList for coord in point]
    localPx = C_compute.WorldTransLocal(flat_paint_joint_M,flatFloatPointList,axis)

    if axis==3:
        weights = C_compute.get_Weights(localPx, px, py, radius, 2)
    else:
        weights = C_compute.get_Weights(localPx,px,py,radius,1)
    skin_cluster.ptw.set(weights)

def get_soft_localPx():
    polygon = get_Geometry(shape_type="mesh")
    skin_cluster = get_skinCluster(polygon)
    mesh = polygon.getShape()

    pm.softSelect(sse=1)
    pm.softSelect(ssc="0,1,1,1,0,1")
    radius = pm.softSelect(q=1, ssd=1)

    pm.softSelect(ssd=radius)
    old_points = mesh.getPoints(space="world")
    flatOldPointsList = [float(coord) for point in old_points for coord in point]

    pm.move([0, 1, 0], r=1)
    new_points = mesh.getPoints(space="world")
    flatNewPointsList = [float(coord) for point in new_points for coord in point]

    pm.move([0, -1, 0], r=1)
    localPx = C_compute.get_Soft_LocalPx(flatNewPointsList,flatOldPointsList)


    return [localPx,skin_cluster]

def paint_soft(px,py,radius):
    templist = get_soft_localPx()
    localPx = templist[0]
    skin_cluster = templist[1]
    weights = C_compute.get_Weights(localPx,px,py,radius,0)
    skin_cluster.ptw.set(weights)

def get_chain_info():
    polygon = get_Geometry(shape_type="mesh")
    skin_cluster = get_skinCluster(polygon)

    # get the order of joints chain
    influences = skin_cluster.getInfluence()
    joints = [joint for joint in influences if not joint.liw.get()]
    unlock_joints = len(joints)
    # create a dict, key is their distance, the value is the pair of joints
    if unlock_joints >1:
        distance_joints = {}
        for i, joint1 in enumerate(joints):
            for joint2 in joints[i + 1:]:
                distance = distance_between_joints(joint1,joint2)
                distance_joints[distance] = [joint1, joint2]
        # find the most close pair of joints
        sorted_joints = distance_joints[min(distance_joints.keys())]
        joints.remove(sorted_joints[0])
        joints.remove(sorted_joints[1])
        # construct the right order of joints chain
        # first create a dict to get the distance between the most close joints and other joints
        for i in range(len(joints)):
            distance ={distance_between_joints(joint1, joint2): [joint1, joint2]
                      for joint1 in joints for joint2 in [sorted_joints[0], sorted_joints[-1]]}
            # find the most close joints pair
            sorted_joint, start_or_end_joint = distance[min(distance.keys())]
            # if star or end joint is the first joint in sorted joints, put it in the start of sorted joints list
            if start_or_end_joint == sorted_joints[0]:
                sorted_joints.insert(0,sorted_joint)
            # or put it in the end of the sorted joints list
            else:
                sorted_joints.append(sorted_joint)
            # then remove the sorted joint from original joint list
            joints.remove(sorted_joint)
        # create an array get the joints index
        joints_index = MIntArray()
        joints_num = len(sorted_joints)
        joints_index.setLength(joints_num)
        for i, joint in enumerate(sorted_joints):
            joints_index[i] = influences.index(joint)
    else:
        sorted_joints = joints
        joints_index = MIntArray()
        joints_index[0] = influences.index(sorted_joints[0])
        joints_num = 1

    # get the original weights before start painting
    points_list = polygon.getShape().getPoints(space="world")
    flat_points_list = [float(coord) for point in points_list for coord in point]
    selected = MSelectionList()
    selected.add(skin_cluster.name())
    selected.add(polygon.getShape().name() + ".vtx[*]")
    print(polygon.getShape().name() + ".vtx[*]")
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
        joint_M = joint.getMatrix(worldSpace=True)
        # convert it to a flat list for computing in C++
        flat_joint_M = [float(value) for row in joint_M for value in row]
        joints_Matrix+=flat_joint_M

    return [flat_points_list,originalWeights_list,joints_Matrix,joints_num,joints_index,components,path,fn_skin]


def paint_Chain_weights(px,py,radius,axis):
    ChaninInfolist = get_chain_info()
    pointlist = ChaninInfolist[0]
    originalWeight_List = ChaninInfolist[1]
    joints_Matrix =ChaninInfolist[2]
    joints_num = ChaninInfolist[3]
    joints_index = ChaninInfolist[4]
    components = ChaninInfolist[5]

    path = ChaninInfolist[6]

    fn_skin = ChaninInfolist[7]


    points_num = int(len(pointlist)/3)
    c_reault_weight = C_compute.get_Chain_Weights(joints_Matrix,pointlist,originalWeight_List,px,py,radius,1,axis)
    api_weights = MDoubleArray()
    api_weights.setLength(len(c_reault_weight))
    for i,w in enumerate(c_reault_weight):
        api_weights[i] = c_reault_weight[i]
    fn_skin.setWeights(path,components,joints_index,api_weights,True)
    refresh_weight()




