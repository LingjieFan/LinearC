gcc -g -c -I ../include/ UnitTest_VectorD.c -o UnitTest_VectorD.o
gcc -g -c -I ../include/ UnitTest_VectorDC.c -o UnitTest_VectorDC.o
gcc -g -c -I ../include/ UnitTest_MatrixD.c -o UnitTest_MatrixD.o
gcc -g -c -I ../include/ UnitTest_MatrixDC.c -o UnitTest_MatrixDC.o
gcc -g -c -I ../include/ UnitTest_LinearC.c -o UnitTest_LinearC.o
gcc -g -c -I ../include/ UnitTest_main.c -o UnitTest_main.o
gcc -g -L ../lib/ UnitTest_main.o UnitTest_LinearC.o UnitTest_MatrixDC.o UnitTest_MatrixD.o UnitTest_VectorD.o UnitTest_VectorDC.o -lObjectC -lLinearC -lblas -llapack -o UnitTest_LinearC.exe
pause
