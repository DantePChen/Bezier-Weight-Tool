#include<Python39/Python/Python.h>
#include"Eigen/Dense"
#include<iostream>

using namespace Eigen;

//list input function
double* pylistOBJ_to_Cpp_Array(PyObject* py_list) {
	int n = PyObject_Length(py_list);
	double* double_array = new double[n];
	for (int i = 0; i < n; i++) {
		double_array[i] = PyFloat_AsDouble(PyList_GetItem(py_list, i));
	}
	return double_array;
}

//list output function
PyObject* Cpp_Array_to_pylistOutput(double* cpp_array, int n) {
	PyObject* pylistResult = PyList_New(n);
	for (int i = 0; i < n; i++) {
		PyList_SetItem(pylistResult, i, PyFloat_FromDouble(cpp_array[i]));
	}
	return pylistResult;
}

double distance(double* jointPostion, double* pointPostion) {
	double distance = pow((pow((jointPostion[0] - pointPostion[0]), 2.0) + pow((jointPostion[1] - pointPostion[1]), 2.0) + pow((jointPostion[2] - pointPostion[2]), 2.0)), 0.5);

	return distance;
}

double bezier(double t, double c_p[4]) {
	return c_p[0] * pow(1 - t, 3.0) + 3.0 * c_p[1] * t * pow(1 - t, 2.0) + 3.0 * c_p[2] * pow(t, 2.0) * (1 - t) + c_p[3] * pow(t, 3.0);
}

static PyObject* py_bezier(PyObject* self, PyObject* args) {
	double c_t;
	double c_p[4];

	if (!PyArg_ParseTuple(args, "d(dddd)",&c_t,&c_p[0], &c_p[1], &c_p[2], &c_p[3])) 
	{
		return NULL;
	}
	double c_result = bezier(c_t, c_p);
	return Py_BuildValue("d", c_result);

}

double bezier_t(double v, double c_p[4]) {
	double min_t = 0;
	double max_t = 1;
	double t;
	double tolerance;
	while (true) {
		t = (min_t + max_t) / 2;
		tolerance = bezier(t, c_p) - v;
		if (tolerance > 0.0001) {
			max_t = t;
		}
		else if (tolerance < -0.0001) {
			min_t = t;
		}
		else {
			return t;
		}

	}
}

static PyObject* py_bezier_t(PyObject* self, PyObject* args) {
	double c_v;
	double c_p[4];

	if (!PyArg_ParseTuple(args, "d(dddd)", &c_v, &c_p[0], &c_p[1], &c_p[2], &c_p[3]))
	{
		return NULL;
	}
	double c_result = bezier(c_v, c_p);
	return Py_BuildValue("d", c_result);

}



double get_pointWeight(double x, double c_px[4], double c_py[4]) {
	double t=0;
	if (x <= 0.0) {
		return c_py[0];
	}
	else if (x >= 1.0) {
		return c_py[3];
	}
	else{
		t = bezier_t(x, c_px);
		return bezier(t,c_py);
	}

}

static PyObject* py_get_pointWeight(PyObject* self, PyObject* args) {
	double c_x;
	double c_px[4];
	double c_py[4];

	if (!PyArg_ParseTuple(args, "d(dddd)(dddd)", &c_x, &c_px[0], &c_px[1], &c_px[2], &c_px[3],
													   &c_py[0], &c_py[1], &c_py[2], &c_py[3]))
	{
		return NULL;
	}
	double c_result = get_pointWeight(c_x, c_px, c_py);
	return Py_BuildValue("d", c_result);

}


double *get_Weights(double* c_localPx, int n, double c_px[4], double c_py[4],double radius,int flag){
	double *c_weights = new double[n];
	if (flag == 1) {
		for (int i = 0; i < n; i++) {
			c_weights[i] = get_pointWeight((c_localPx[i] / radius/2) + 0.5, c_px, c_py);
		}
	}
	else if (flag == 2) {
		for (int i = 0; i < n; i++) {
			c_weights[i] = get_pointWeight(1-(c_localPx[i] /radius), c_px, c_py);
		}
	}
	else{
		for (int i = 0; i < n; i++) {
			c_weights[i] = get_pointWeight(c_localPx[i], c_px, c_py);
		}
	}
	
	return c_weights;
}



// deal with the python inputs and outputs
static PyObject* py_get_Weights(PyObject* self, PyObject* args) {
	// declare C++ args
	double c_px[4];
	double c_py[4];
	double radius;
	int flag;

	// convert from Python args to C++ args
	PyObject *py_localPx;
	if (!PyArg_ParseTuple(args, "O(dddd)(dddd)di", &py_localPx,
		&c_px[0], &c_px[1], &c_px[2], &c_px[3],
		&c_py[0], &c_py[1], &c_py[2], &c_py[3],&radius,&flag)
		) {
		return NULL;
	}
		
	int n = PyObject_Length(py_localPx);
	// create C++ array from Python list
	// declare a C++ array and write the value from Python list
	double* c_localPx = pylistOBJ_to_Cpp_Array(py_localPx);

	// get result form C++ function,this is a array
	double* c_weights = get_Weights(c_localPx, n, c_px, c_py, radius,flag);

	// convert C++ array to Python list
	// declare a python list
	// write the value from C++ array
	PyObject* py_weights = Cpp_Array_to_pylistOutput(c_weights, n);

	//free memory
	delete[] c_localPx;
	delete[] c_weights;

	return py_weights;

}

double* WorldTransLocal(double* c_jointMatrix, double* c_pointList, int n,int axis) {
	int count = n / 3;
	double* resultArray = new double[count];
	if (axis == 3) {
		double* jointPostion = new double[3];
		for (int i = 0; i < 3; i++) {
			jointPostion[i] = c_jointMatrix[12+i];
		}
		for (int i = 0; i < count; i++) {
			double* pointPostion = new double[3];
			for (int j = 0; j < 3; j++) {
				pointPostion[j] = c_pointList[3 * i + j];
			}

			resultArray[i] = distance(jointPostion, pointPostion);

		}
		// find the min of whole array
		int min = resultArray[0];
		for (int i = 1; i < count; i++) {
			if (resultArray[i] < min) {
				min = resultArray[i];
			}
		}
		for (int i = 0; i < count; i++) {
			resultArray[i] -= min;
		}


	}
	else {
		//convert C++ array to Eigen Matrix
		MatrixXd jointMatrix = MatrixXd::Zero(4, 4);
		MatrixXd inverseJointMatrix = MatrixXd::Zero(4, 4);
		for (int i = 0; i < 4; i++) {
			jointMatrix(i, 0) = c_jointMatrix[4 * i];
			jointMatrix(i, 1) = c_jointMatrix[4 * i + 1];
			jointMatrix(i, 2) = c_jointMatrix[4 * i + 2];
			jointMatrix(i, 3) = c_jointMatrix[4 * i + 3];
		}

		inverseJointMatrix = jointMatrix.inverse();
		for (int i = 0; i < count; i++) {
			MatrixXd pointMatrix = MatrixXd::Zero(1, 4);
			pointMatrix(0, 0) = c_pointList[3 * i];
			pointMatrix(0, 1) = c_pointList[3 * i + 1];
			pointMatrix(0, 2) = c_pointList[3 * i + 2];
			pointMatrix(0, 3) = 1;
			MatrixXd resultMatrix = MatrixXd::Zero(1, 4);
			resultMatrix = pointMatrix * inverseJointMatrix;
			if (axis == 0) {
				resultArray[i] = resultMatrix(0, 0);
			}
			else if (axis == 1) {
				resultArray[i] = resultMatrix(0, 1);
			}
			else {
				resultArray[i] = resultMatrix(0, 2);
			}

		}

	}
	

	return resultArray;
	
}

static PyObject* py_WorldTransLocal(PyObject* self, PyObject* args) {
	PyObject *py_jointMatrix;
	PyObject *py_pointList;
	int axis;
	if (!PyArg_ParseTuple(args, "OOi", &py_jointMatrix, &py_pointList,&axis)) 
	{
		return NULL;
	}

	double* c_jointMatrix = pylistOBJ_to_Cpp_Array(py_jointMatrix);


	int n2 = PyObject_Length(py_pointList);
	double* c_pointList = pylistOBJ_to_Cpp_Array(py_pointList);
	

	double* c_local = WorldTransLocal(c_jointMatrix, c_pointList,n2,axis);
	int count = n2 / 3;
	// write the value from C++ array
	PyObject* py_local = Cpp_Array_to_pylistOutput(c_local, count);
	
	//free memory
	delete[] c_jointMatrix;
	delete[] c_pointList;
	delete[] c_local;

	return py_local;
}


double* get_Soft_LocalPx(double* c_new_points, double* c_old_points,int n) {
	double* c_localPx = new double[n];
	for (int i = 0; i < n; i++) {
		c_localPx[i] = (c_new_points[3*i+1] - c_old_points[3*i+1]);

	}
	return c_localPx;
}

static PyObject* py_get_Soft_LocalPx(PyObject* self, PyObject* args) {
	PyObject* py_new_points;
	PyObject* py_old_points;
	if (!PyArg_ParseTuple(args, "OO", &py_new_points, &py_old_points))
	{
		return NULL;
	}
	int n = PyObject_Length(py_new_points)/3;
	double* c_new_points = pylistOBJ_to_Cpp_Array(py_new_points);
	double* c_old_points = pylistOBJ_to_Cpp_Array(py_old_points);

	double* c_LocalPx = get_Soft_LocalPx(c_new_points, c_old_points, n);
	
	PyObject* py_localPx = Cpp_Array_to_pylistOutput(c_LocalPx, n);
	// free memory
	delete[] c_new_points;
	delete[] c_old_points;
	delete[] c_LocalPx;

	return py_localPx;

}

double* get_original_weight(double* c_orignalWeightList, int joint_num, int point_num) {
	
	
	double* totalResult = new double[point_num];
	for (int j = 0; j < point_num; j++) {
		for (int i = 0; i < joint_num; i++) {
			totalResult[j] += c_orignalWeightList[j* joint_num +i];
		}

	}

	for (int i = 0; i < point_num; i++) {
		if (totalResult[i] >= 0.9999) {
			totalResult[i] = 1.0;
		}
	}

	return totalResult;
}

double* get_Chain_Weights(double* jointsMatrix, double* pointlist,double* originalWeight, 
					  double c_px[4], double c_py[4],double radius, const int joint_num,
	                 const int points_num,int flag, int axis) {
	double* raw_weights = new double[joint_num * points_num];
	double* result_weights = new double[joint_num * points_num];
	int n = points_num * 3;
	
	// caculate weights for every joint
	for (int k = 0; k < joint_num; k++) {
		//seperate every the world matrix of every joint
		double* tempJointMatrixArray = new double[16];
		for (int i = 0; i < 16; i++) {
			tempJointMatrixArray[i] = jointsMatrix[k * 16 + i];
		}

		// calculate the except weights for every joints
		double* localPx = WorldTransLocal(tempJointMatrixArray, pointlist, n,axis);
		double* weight = get_Weights(localPx, points_num, c_px, c_py, radius, flag);

		// Copy the calculated weights to the raw_weights array
		for (int i = 0; i < points_num; i++) {
			raw_weights[k * points_num + i] = weight[i];
		}
		delete[] tempJointMatrixArray;
		delete[] localPx;
		delete[] weight;
	}

	for (int i = 0; i < points_num; i++) {
		raw_weights[i] = 1.0;

	}

	// write the value to the result_weights
	for (int i = 0; i < points_num * joint_num; i++) {
		result_weights[i] = raw_weights[i];

	}
	
	for (int j = 0; j < joint_num - 1; j++) {
		for (int i = 0; i < points_num; i++) {
			int joint1 = j * points_num + i;
			int joint2 = (j + 1) * points_num + i;
			result_weights[i*joint_num+j] = raw_weights[joint1] - raw_weights[joint2];
		}

	}

	for (int i = 0; i < points_num; i++) {
		result_weights[i * joint_num + joint_num - 1] = raw_weights[(joint_num - 1) * points_num + i];
	}
	

	for (int j = 0; j < points_num; j++) {
		for (int i = 0; i < joint_num; i++) {
			result_weights[i * joint_num + j] *= originalWeight[j];
		}

	}

	return result_weights;
}

static PyObject* py_get_Chain_Weights(PyObject* self, PyObject* args) {
	double c_px[4];
	double c_py[4];
	double radius;
	int flag;
	int axis;

	PyObject* py_jointsMatrix;
	PyObject* py_pointlist;
	PyObject* py_originalWeightList;
	if (!PyArg_ParseTuple(args, "OOO(dddd)(dddd)dii", &py_jointsMatrix, &py_pointlist,&py_originalWeightList,
		& c_px[0], &c_px[1], &c_px[2], &c_px[3],
		&c_py[0], &c_py[1], &c_py[2], &c_py[3], &radius, &flag,&axis))
	{
		return NULL;
	}

	const int joint_num = PyObject_Length(py_jointsMatrix) / 16;
	int joint_matrix = PyObject_Length(py_jointsMatrix);
	const int point_num = PyObject_Length(py_pointlist)/3;
	int total_elements = PyObject_Length(py_originalWeightList);
	double* c_jointsMatrix = pylistOBJ_to_Cpp_Array(py_jointsMatrix);
	double* c_pointlist = pylistOBJ_to_Cpp_Array(py_pointlist);
	double* c_originalWeightList = pylistOBJ_to_Cpp_Array(py_originalWeightList);

	double* originalWeights = get_original_weight(c_originalWeightList, joint_num, point_num);
	double* result_weights = get_Chain_Weights(c_jointsMatrix, c_pointlist, originalWeights, c_px, c_py, radius, joint_num, point_num, flag, axis);

	PyObject* py_result_weights = Cpp_Array_to_pylistOutput(result_weights, total_elements);

	//free memory
	delete[] c_jointsMatrix;
	delete[] c_pointlist;
	delete[] c_originalWeightList;
	delete[] originalWeights;
	delete[] result_weights;

	return py_result_weights;
	

}



static PyMethodDef BezierWeightToolComputeMethods[] = {
	{"get_Weights", py_get_Weights, METH_VARARGS, "Execute a shell command"},
	{"bezier", py_bezier, METH_VARARGS, "Execute a shell command"},
	{"bezier_t", py_bezier_t, METH_VARARGS, "Execute a shell command"},
	{"get_pointWeight", py_get_pointWeight, METH_VARARGS, "Execute a shell command"},
	{"WorldTransLocal", py_WorldTransLocal, METH_VARARGS, "Execute a shell command"},
	{"get_Soft_LocalPx", py_get_Soft_LocalPx, METH_VARARGS, "Execute a shell command"},
	{"get_Chain_Weights", py_get_Chain_Weights, METH_VARARGS, "Execute a shell command"},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef BezierWeightToolCompute = {
	PyModuleDef_HEAD_INIT,
	"BezierWeightToolCompute", // name of module
	NULL,                      // docstring (set to NULL if not used)
	-1,
	BezierWeightToolComputeMethods
};

PyMODINIT_FUNC PyInit_BezierWeightToolCompute(void) {
	return PyModule_Create(&BezierWeightToolCompute);
}
