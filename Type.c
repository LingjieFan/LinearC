#include "Type.h"

static Type _OBJECT = {NULL};

const Type *OBJECT = & _OBJECT;

static Type _TENSOR = {&_OBJECT};

const Type *TENSOR = & _TENSOR;

static Type _NUM = {&_TENSOR};

static Type _NUMD = {&_NUM};

static Type _NUMDC = {&_NUM};

const Type *NUM = & _NUM;

const Type *NUMD = & _NUMD;

const Type *NUMDC = & _NUMDC;

static Type _VECTOR = {&_TENSOR};

static Type _VECTORD = {&_VECTOR};

static Type _VECTORDC = {&_VECTOR};

const Type *VECTOR = & _VECTOR;

const Type *VECTORD = & _VECTORD;

const Type *VECTORDC = & _VECTORDC;

static Type _MATRIX = {&_TENSOR};

static Type _MATRIXD = {&_MATRIX};

static Type _MATRIXDC = {&_MATRIX};

const Type *MATRIX = & _MATRIX;

const Type *MATRIXD = & _MATRIXD;

const Type *MATRIXDC = & _MATRIXDC;
