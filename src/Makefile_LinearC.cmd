gcc -g -c -fPIC -I ../include MatrixD.c -o MatrixD.o
gcc -g -c -fPIC -I ../include MatrixDC.c -o MatrixDC.o
gcc -g -c -fPIC -I ../include VectorD.c -o VectorD.o
gcc -g -c -fPIC -I ../include VectorDC.c -o VectorDC.o
gcc -g -shared -fPIC -L ../lib MatrixD.o MatrixDC.o VectorD.o VectorDC.o -lObjectC -lblas -llapack -o ../lib/LinearC.dll
pause
