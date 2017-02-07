#include <python2.7/Python.h>
#include <string>

void test() {
  Py_Initialize();
  if (!Py_IsInitialized()) {
    fprintf(stderr, "can't initialize python environment.\n");
    exit(1);
  }

  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('/home/zengping/work/code/waterbox/QMCPP/QM/python')");

  PyObject* pModule = PyImport_ImportModule("test_calculations_wtr_new");
  PyObject* pDict = PyModule_GetDict(pModule);
  if(!pDict) {
    fprintf(stderr, "can't find dictionary.\n");
    exit(1);
  }

  PyObject* pClassQM = PyDict_GetItemString(pDict, "QMInterpolation");
  if(!pClassQM) {
    fprintf(stderr, "can't find QMInterpolation class.\n");
    exit(1);
  }

  // Construct QM Instance
  PyObject* tmp = PyString_FromString("wtr");
  PyObject* pInstanceArgs = PyTuple_Pack(1, tmp);
  if (!pInstanceArgs) {
    PyErr_Print();
    exit(1);
  }
  Py_DECREF(tmp);

  PyObject* pInstanceQM = PyInstance_New(pClassQM, pInstanceArgs, NULL);
  if(!pInstanceQM)
  {
    PyErr_Print();
    fprintf(stderr, "can't create QM instance.\n");
    exit(1);
  }
  PyObject* lhs = Py_BuildValue("[(s, [dddd]), (s, [dddd]), (s, [dddd])]",
                                "O", 8.0, 0.10293690, 0.11093940, 0.07435482,
                                "H", 1.0, 0.73097151, 0.83572755, 0.02040124,
                                "H", 1.0, -0.22532655, 0.13689786, 0.97669930);
  if(!lhs) {
    PyErr_Print();
    exit(1);
  }

  PyObject* rhs = Py_BuildValue("[(s, [dddd]), (s, [dddd]), (s, [dddd])]",
                                "O", 8.0, 2.85534030, 0.60573860, -1.03993027,
                                "H", 1.0, 1.93176673, 0.51127346, -0.79346626,
                                "H", 1.0, 3.16776562, -0.29302760, -1.17133059);
  if(!rhs) {
    PyErr_Print();
    exit(1);
  }
  PyObject *pReturn = NULL;
  pReturn = PyObject_CallMethod(pInstanceQM,
                                (char*)"calculate",
                                (char*)"(OO)", lhs, rhs);

  if (!pReturn) {
    PyErr_Print();
    exit(1);
  }
  double tempe;
  double r;
  if (!PyArg_ParseTuple(pReturn, "dd", &tempe, &r)) {
    fprintf(stderr, "Failed to parse the return val\n");
    exit(1);
  }

  printf("tempe:%lf  r: %lf\n", tempe, r);

  Py_DECREF(lhs);
  Py_DECREF(rhs);
  Py_DECREF(pReturn);
  Py_DECREF(pInstanceArgs);
  Py_DECREF(pInstanceQM);
  Py_DECREF(pClassQM);
  Py_DECREF(pDict);
  Py_DECREF(pModule);
  Py_Finalize();
}

int main(void) {
  test();
  return 0;
}
